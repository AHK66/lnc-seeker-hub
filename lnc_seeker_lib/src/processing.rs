// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::collections::{HashMap, HashSet, BTreeSet};

use crate::models::ReadInfo;

pub fn get_junction_reads(
    bam_path: &str,
    reference: &str,
    start_target: usize,
    end_target: usize,
    min_mq: u8,
    max_reads: usize,
    _genome_path: &Option<String>,
    filter_clean: bool,
) -> Result<Vec<ReadInfo>, Box<dyn std::error::Error + Send + Sync>> {
    let read_provider = crate::reads_manager::get_read_provider();
    
    // Ensure cache is loaded from disk if available but not in RAM
    if !read_provider.is_in_ram(bam_path) {
        let _ = read_provider.load_from_disk(bam_path);
    }

    if let Some((filtered_reads, status)) = read_provider.get_filtered_reads(bam_path, reference, start_target, end_target, min_mq) {
        if status == crate::reads_manager::CacheStatus::Complete {
            println!("Rust: Using {} cached reads for {}", filtered_reads.len(), bam_path);
            let mut name_to_read: HashMap<String, ReadInfo> = HashMap::new();
            let mut potential_mates: HashMap<String, Vec<crate::models::ReadSegment>> = HashMap::new();
            let mut keep_names = HashSet::new();

            for info in filtered_reads {
                let mut has_target_junction = false;
                for i in 0..info.segments.len().saturating_sub(1) {
                    if !info.segments[i].is_followed_by_deletion && info.segments[i].end == start_target && info.segments[i+1].start == end_target {
                        if filter_clean {
                            if info.segments[i].mismatches == 0 && info.segments[i].insertions == 0 && 
                               info.segments[i+1].mismatches == 0 && info.segments[i+1].insertions == 0 {
                                has_target_junction = true;
                                break;
                            }
                        } else {
                            has_target_junction = true;
                            break;
                        }
                    }
                }

                if has_target_junction {
                    keep_names.insert(info.name.clone());
                    let target = name_to_read.entry(info.name.clone()).or_insert(ReadInfo {
                        name: info.name.clone(),
                        reference: info.reference.clone(),
                        mapping_quality: info.mapping_quality,
                        strand: info.strand,
                        start: info.start,
                        end: info.end,
                        segments: Vec::new(),
                    });
                    target.segments.extend(info.segments.iter().cloned());
                    target.start = target.start.min(info.start);
                    target.end = target.end.max(info.end);

                    if let Some(mate_segs) = potential_mates.remove(&info.name) {
                        for s in &mate_segs {
                            target.start = target.start.min(s.start);
                            target.end = target.end.max(s.end);
                        }
                        target.segments.extend(mate_segs);
                    }
                } else if keep_names.contains(&info.name) {
                    if let Some(target) = name_to_read.get_mut(&info.name) {
                        for s in &info.segments {
                            target.start = target.start.min(s.start);
                            target.end = target.end.max(s.end);
                        }
                        target.segments.extend(info.segments.iter().cloned());
                    }
                } else {
                    // Check for bridge
                    let mut mate_anchored_start = false;
                    let mut mate_anchored_end = false;
                    for s in &info.segments {
                        if s.end == start_target { 
                            if !filter_clean || (s.mismatches == 0 && s.insertions == 0) {
                                mate_anchored_start = true; 
                            }
                        }
                        if s.start == end_target { 
                            if !filter_clean || (s.mismatches == 0 && s.insertions == 0) {
                                mate_anchored_end = true; 
                            }
                        }
                    }

                    if let Some(existing_segs) = potential_mates.get(&info.name) {
                        for s in existing_segs {
                            if s.end == start_target { 
                                if !filter_clean || (s.mismatches == 0 && s.insertions == 0) {
                                    mate_anchored_start = true; 
                                }
                            }
                            if s.start == end_target { 
                                if !filter_clean || (s.mismatches == 0 && s.insertions == 0) {
                                    mate_anchored_end = true; 
                                }
                            }
                        }
                    }

                    if mate_anchored_start && mate_anchored_end {
                        keep_names.insert(info.name.clone());
                        let mut all_segs = info.segments.clone();
                        if let Some(extracted) = potential_mates.remove(&info.name) {
                            all_segs.extend(extracted);
                        }
                        let min_s = all_segs.iter().map(|s| s.start).min().unwrap_or(info.start);
                        let max_e = all_segs.iter().map(|s| s.end).max().unwrap_or(info.end);

                        name_to_read.insert(info.name.clone(), ReadInfo {
                            name: info.name.clone(),
                            reference: info.reference.clone(),
                            mapping_quality: info.mapping_quality,
                            strand: info.strand,
                            start: min_s,
                            end: max_e,
                            segments: all_segs,
                        });
                    } else {
                        potential_mates.entry(info.name.clone()).or_default().extend(info.segments.iter().cloned());
                    }
                }
            }
            return Ok(finalize_results(name_to_read, max_reads));
        }
    }

    println!("Rust: Cache not available for {}, and BAM fallback is disabled.", bam_path);
    Err(format!("Compact storage cache (.lnc_cache.bin) missing or incomplete for {}. Please re-run lnc-seeker-collect.", bam_path).into())
}

fn finalize_results(name_to_read: HashMap<String, ReadInfo>, max_reads: usize) -> Vec<ReadInfo> {
    // Sort by largest span between outermost coordinates descending
    let mut results: Vec<ReadInfo> = name_to_read.into_values().collect();
    results.sort_by(|a, b| {
        let span_a = a.end.saturating_sub(a.start);
        let span_b = b.end.saturating_sub(b.start);
        span_b.cmp(&span_a)
    });
    
    // Take top N
    if results.len() > max_reads {
        results.truncate(max_reads);
    }

    // Sort segments for each read to ensure UI draws them correctly
    for info in &mut results {
        info.segments.sort_by_key(|s| s.start);
    }
    
    results
}


pub fn downsample_coverage(
    positions: Vec<usize>,
    depths_bg: Vec<f64>,
    depths_fg: Vec<f64>,
    target_points: usize,
    forced_positions: Vec<usize>,
) -> (Vec<usize>, Vec<f64>, Vec<f64>) {
    let n = positions.len();
    if n <= target_points || target_points < 10 {
        return (positions, depths_bg, depths_fg);
    }

    let mut selected_indices = BTreeSet::new();
    selected_indices.insert(0);
    selected_indices.insert(n - 1);

    // Force include sharp cliffs identified as junction points
    for &pos in &forced_positions {
        if let Ok(idx) = positions.binary_search(&pos) {
            selected_indices.insert(idx);
            // Include neighbor to preserve the vertical edge of the cliff
            if idx > 0 {
                selected_indices.insert(idx - 1);
            }
            if idx + 1 < n {
                selected_indices.insert(idx + 1);
            }
        }
    }

    // We want to achieve roughly target_points output.
    let num_bins = (target_points / 6).max(1);
    let bin_size = n as f64 / num_bins as f64;

    for i in 0..num_bins {
        let start = (i as f64 * bin_size) as usize;
        let end = (((i + 1) as f64 * bin_size) as usize).min(n);
        if start >= end {
            continue;
        }

        let mut min_bg_idx = start;
        let mut max_bg_idx = start;
        let mut min_fg_idx = start;
        let mut max_fg_idx = start;
        let mut max_diff_idx = start;
        let mut max_diff = -1.0;

        for j in start..end {
            if depths_bg[j] < depths_bg[min_bg_idx] {
                min_bg_idx = j;
            }
            if depths_bg[j] > depths_bg[max_bg_idx] {
                max_bg_idx = j;
            }
            if depths_fg[j] < depths_fg[min_fg_idx] {
                min_fg_idx = j;
            }
            if depths_fg[j] > depths_fg[max_fg_idx] {
                max_fg_idx = j;
            }

            if j > 0 {
                let diff_bg = (depths_bg[j] - depths_bg[j - 1]).abs();
                let diff_fg = (depths_fg[j] - depths_fg[j - 1]).abs();
                let diff = diff_bg.max(diff_fg);
                if diff > max_diff {
                    max_diff = diff;
                    max_diff_idx = j;
                }
            }
        }

        selected_indices.insert(min_bg_idx);
        selected_indices.insert(max_bg_idx);
        selected_indices.insert(min_fg_idx);
        selected_indices.insert(max_fg_idx);
        selected_indices.insert(max_diff_idx);
        if max_diff_idx > 0 {
            selected_indices.insert(max_diff_idx - 1);
        }
    }

    let mut res_x = Vec::with_capacity(selected_indices.len());
    let mut res_bg = Vec::with_capacity(selected_indices.len());
    let mut res_fg = Vec::with_capacity(selected_indices.len());

    for idx in selected_indices {
        res_x.push(positions[idx]);
        res_bg.push(depths_bg[idx]);
        res_fg.push(depths_fg[idx]);
    }

    (res_x, res_bg, res_fg)
}

pub fn get_junction_reads_batch(
    bam_path: &str,
    reference: &str,
    junctions: &[(usize, usize)],
    min_mq: u8,
    max_reads: usize,
    _genome_path: &Option<String>,
    filter_clean: bool,
) -> Result<HashMap<String, Vec<ReadInfo>>, Box<dyn std::error::Error + Send + Sync>> {
    let read_provider = crate::reads_manager::get_read_provider();

    // Ensure cache is loaded from disk if available but not in RAM
    if !read_provider.is_in_ram(bam_path) {
        let _ = read_provider.load_from_disk(bam_path);
    }

    let mut results_map = HashMap::new();

    if let Some((batch_filtered, _status)) =
        read_provider.get_filtered_reads_batch(bam_path, reference, junctions, min_mq)
    {
        for (&(js, je), filtered_reads) in &batch_filtered {
            let mut name_to_read: HashMap<String, ReadInfo> = HashMap::new();
            let mut potential_mates: HashMap<String, Vec<crate::models::ReadSegment>> = HashMap::new();
            let mut keep_names = HashSet::new();

            for info in filtered_reads {
                let mut has_target_junction = false;
                for i in 0..info.segments.len().saturating_sub(1) {
                    if !info.segments[i].is_followed_by_deletion
                        && info.segments[i].end == js
                        && info.segments[i + 1].start == je
                    {
                        if filter_clean {
                            if info.segments[i].mismatches == 0
                                && info.segments[i].insertions == 0
                                && info.segments[i + 1].mismatches == 0
                                && info.segments[i + 1].insertions == 0
                            {
                                has_target_junction = true;
                                break;
                            }
                        } else {
                            has_target_junction = true;
                            break;
                        }
                    }
                }

                if has_target_junction {
                    keep_names.insert(info.name.clone());
                    let target = name_to_read.entry(info.name.clone()).or_insert(ReadInfo {
                        name: info.name.clone(),
                        reference: info.reference.clone(),
                        mapping_quality: info.mapping_quality,
                        strand: info.strand,
                        start: info.start,
                        end: info.end,
                        segments: Vec::new(),
                    });
                    target.segments.extend(info.segments.iter().cloned());
                    target.start = target.start.min(info.start);
                    target.end = target.end.max(info.end);

                    if let Some(mate_segs) = potential_mates.remove(&info.name) {
                        for s in &mate_segs {
                            target.start = target.start.min(s.start);
                            target.end = target.end.max(s.end);
                        }
                        target.segments.extend(mate_segs);
                    }
                } else if keep_names.contains(&info.name) {
                    if let Some(target) = name_to_read.get_mut(&info.name) {
                        for s in &info.segments {
                            target.start = target.start.min(s.start);
                            target.end = target.end.max(s.end);
                        }
                        target.segments.extend(info.segments.iter().cloned());
                    }
                } else {
                    let mut mate_anchored_start = false;
                    let mut mate_anchored_end = false;
                    for s in &info.segments {
                        if s.end == js {
                            if !filter_clean || (s.mismatches == 0 && s.insertions == 0) {
                                mate_anchored_start = true;
                            }
                        }
                        if s.start == je {
                            if !filter_clean || (s.mismatches == 0 && s.insertions == 0) {
                                mate_anchored_end = true;
                            }
                        }
                    }

                    if let Some(existing_segs) = potential_mates.get(&info.name) {
                        for s in existing_segs {
                            if s.end == js {
                                if !filter_clean || (s.mismatches == 0 && s.insertions == 0) {
                                    mate_anchored_start = true;
                                }
                            }
                            if s.start == je {
                                if !filter_clean || (s.mismatches == 0 && s.insertions == 0) {
                                    mate_anchored_end = true;
                                }
                            }
                        }
                    }

                    if mate_anchored_start && mate_anchored_end {
                        keep_names.insert(info.name.clone());
                        let target = name_to_read.entry(info.name.clone()).or_insert(ReadInfo {
                            name: info.name.clone(),
                            reference: info.reference.clone(),
                            mapping_quality: info.mapping_quality,
                            strand: info.strand,
                            start: info.start,
                            end: info.end,
                            segments: Vec::new(),
                        });
                        target.segments.extend(info.segments.iter().cloned());
                        target.start = target.start.min(info.start);
                        target.end = target.end.max(info.end);

                        if let Some(mate_segs) = potential_mates.remove(&info.name) {
                            for s in &mate_segs {
                                target.start = target.start.min(s.start);
                                target.end = target.end.max(s.end);
                            }
                            target.segments.extend(mate_segs);
                        }
                    } else {
                        potential_mates
                            .entry(info.name.clone())
                            .or_insert_with(Vec::new)
                            .extend(info.segments.iter().cloned());
                    }
                }
            }

            let mut final_reads: Vec<ReadInfo> = name_to_read.into_values().collect();
            for r in &mut final_reads {
                r.segments.sort_by_key(|s| s.start);
            }
            final_reads.sort_by_key(|r| r.start);
            if final_reads.len() > max_reads {
                final_reads.truncate(max_reads);
            }

            results_map.insert(format!("{}-{}", js, je), final_reads);
        }
    }

    Ok(results_map)
}
