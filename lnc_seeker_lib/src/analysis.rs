// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader};
use std::time::Instant;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use noodles::{bam, sam};
use sam::alignment::Record;
use rayon::prelude::*;

use crate::config::Config;
use crate::models::{AnalysisResult, SampleResult};
use crate::progress::ProgressData;
use crate::coverage::{Coverage, JunctionStore};
use crate::assembly::AssemblyReport;
use crate::annotations::get_annotations;

pub fn cleanup_csv_files() -> io::Result<()> {
    let patterns = [
        "_junction_points.csv",
        "_junction_spans.csv",
        "_junction_points_hq.csv",
        "_junction_spans_hq.csv",
        "_coverage.csv",
    ];

    let entries = std::fs::read_dir(".")?;
    for entry in entries {
        let entry = entry?;
        let path = entry.path();
        if path.is_file() {
            if let Some(file_name) = path.file_name().and_then(|s| s.to_str()) {
                let should_delete = file_name == "annotations.csv" || patterns.iter().any(|&p| file_name.ends_with(p));
                if should_delete {
                    let _ = std::fs::remove_file(path);
                }
            }
        }
    }
    Ok(())
}

pub fn load_assembly_report(path: &str) -> io::Result<AssemblyReport> {
    AssemblyReport::parse_from_file(path)
}

pub fn run_analysis(config: &Config, progress: Arc<ProgressData>) -> io::Result<AnalysisResult> {
    let t_total = Instant::now();
    let res = run_analysis_inner(config, Arc::clone(&progress));
    
    let read_provider = crate::reads_manager::get_read_provider();
    read_provider.benchmark.record_stage("Total Pipeline", t_total.elapsed());
    read_provider.benchmark.print_report();
    let _ = read_provider.benchmark.write_csv("benchmarking/data/benchmark_analysis.csv");

    progress.stage.store(0, Ordering::SeqCst);
    res
}

pub fn populate_cache_from_records<R>(
    bam_path: &str,
    records: &[R],
    header: &sam::Header,
    genome_provider: &mut crate::genome::GenomeProvider,
    read_provider: &crate::reads_manager::ReadProvider,
) -> io::Result<()> 
where R: Record {
    let mut read_infos = Vec::with_capacity(records.len());
    let mut ref_names = Vec::new();
    let mut ref_lengths = Vec::new();
    
    for (name, seq) in header.reference_sequences().iter() {
        ref_names.push(String::from_utf8_lossy(name.as_ref()).to_string());
        ref_lengths.push(seq.length().get());
    }

    for record in records {
        let ref_id = if let Some(id) = record.reference_sequence_id(header).transpose()? {
             Some(id)
        } else {
             None
        };

        let ref_seq = if let Some(ref_id) = ref_id {
            if let Some((ref_name_b, _)) = header.reference_sequences().get_index(ref_id) {
                let ref_name = std::str::from_utf8(ref_name_b.as_ref()).unwrap_or("");
                genome_provider.get_sequence(ref_name).unwrap_or(None).map(|s| s.as_ref())
            } else {
                None
            }
        } else {
            None
        };

        if let Some(info) = crate::models::ReadInfo::from_record(record, header, ref_seq)? {
            read_infos.push(info);
        }
    }

    read_provider.commit_caching_compact(
        bam_path,
        read_infos,
        ref_names,
        ref_lengths,
        crate::reads_manager::CacheStatus::Complete,
    ).map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;

    Ok(())
}

pub fn populate_bam_cache(
    src: &str,
    progress: &Arc<ProgressData>,
    read_provider: &crate::reads_manager::ReadProvider,
    genome_path: &Option<String>,
) -> io::Result<()> {
    let t_start = Instant::now();
    println!("Rust: Populating cache from BAM: {}", src);
    
    let file = File::open(src)?;
    let mut reader = bam::io::Reader::new(noodles::bgzf::io::Reader::new(BufReader::new(file)));
    let header = reader.read_header()?;
    
    let mut genome_provider = crate::genome::GenomeProvider::new(genome_path)?;
    
    let mut ref_names = Vec::new();
    let mut ref_lengths = Vec::new();
    for (name, seq) in header.reference_sequences().iter() {
        ref_names.push(name.to_string());
        ref_lengths.push(seq.length().get());
    }

    let mut builder = crate::reads_manager::IncrementalCacheBuilder::new(src, ref_names, ref_lengths);

    // Incremental memory estimation state
    let mut est_bytes: usize = 0;
    let limit_mb = *read_provider.max_memory_mb.read().unwrap_or_else(|e| e.into_inner());
    let limit_bytes = (limit_mb * 1024.0 * 1024.0) as usize;
    let current_usage_bytes = (read_provider.get_total_memory_usage_mb() * 1024.0 * 1024.0) as usize;

    let start_loop = Instant::now();
    let mut count = 0;
    for result in reader.records() {
        let record = result?;

        let ref_id_res = record.reference_sequence_id();
        let ref_id = if let Some(id) = ref_id_res.transpose()? {
             Some(id)
        } else {
             None
        };

        let ref_seq = if let Some(ref_id) = ref_id {
            if let Some((ref_name_b, _)) = header.reference_sequences().get_index(ref_id) {
                let ref_name = std::str::from_utf8(ref_name_b.as_ref()).unwrap_or("");
                genome_provider.get_sequence(ref_name).unwrap_or(None).map(|s| s.as_ref())
            } else {
                None
            }
        } else {
            None
        };

        if let Some(info) = crate::models::ReadInfo::from_record(&record, &header, ref_seq)? {
            // Incremental estimation: CompactRead(24) + segments(8) + name(48 overhead)
            let read_overhead = 24 + (info.segments.len() * 8) + (info.name.len() + 48);
            est_bytes += read_overhead;

            // Check every 1000 records
            if count % 1000 == 0 && current_usage_bytes + est_bytes > limit_bytes {
                let msg = format!("Rust: Global cache memory limit reached ({:.1} MB + est {:.1} MB > {:.1} MB). Population aborted for {}.", 
                            current_usage_bytes as f64 / 1024.0 / 1024.0, 
                            est_bytes as f64 / 1024.0 / 1024.0, 
                            limit_mb, src);
                return Err(io::Error::new(io::ErrorKind::Other, msg));
            }

            builder.add_read(info);
        }

        count += 1;
        if count % 10000 == 0 {
            progress.current.fetch_add(10000, Ordering::Relaxed);
        }
    }
    
    println!("[BENCHMARK] Rust: Record population ({} records, ~{:.1} MB) took {:?}", 
             count, est_bytes as f64 / 1024.0 / 1024.0, start_loop.elapsed());
    
    read_provider.commit_builder(
        src, 
        builder, 
        crate::reads_manager::CacheStatus::Complete
    )?;

    read_provider.benchmark.record_stage(&format!("Populate Cache: {}", src), t_start.elapsed());

    Ok(())
}

fn run_analysis_inner(config: &Config, progress: Arc<ProgressData>) -> io::Result<AnalysisResult> {
    let t_preamble = Instant::now();
    let mut all_covered_ranges: HashMap<String, (usize, usize)> = HashMap::new();
    let bam_paths = config.get_paths();
    println!("********** Rust: Starting analysis for {} BAM files.", bam_paths.len());

    let read_provider = crate::reads_manager::get_read_provider();
    read_provider.benchmark.observe_peak();

    // Purge caches for samples that are NO LONGER in the current session
    read_provider.retain_only_selected(&bam_paths);

    // Update global cache limit if provided in config
    if let Some(limit_mb) = config.data_selection.max_cache_memory_mb {
        read_provider.set_max_memory_mb(limit_mb);
        println!("Rust: Set global cache memory limit to {:.1} MB", limit_mb);
    }
    
    // Reset progress
    progress.stage.store(1, Ordering::SeqCst); // 1: BAM Analysis
    progress.current.store(0, Ordering::SeqCst);
    progress.total.store(0, Ordering::SeqCst);

    read_provider.benchmark.record_stage("BAM Analysis Preamble", t_preamble.elapsed());

    let results: Vec<io::Result<(String, SampleResult, HashMap<String, (usize, usize)>)>> = bam_paths.par_iter().map(|src| {
        let progress = Arc::clone(&progress);
        let start_time = Instant::now();

        let read_provider = crate::reads_manager::get_read_provider();
        
        let sample_name = config.data_selection.bam_to_cohort.get(src)
            .cloned()
            .unwrap_or_else(|| {
                std::path::Path::new(src)
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("sample")
                    .to_string()
            });

        read_provider.benchmark.observe_cohort_peak(&sample_name);

        let mut junctions = JunctionStore::default();
        let mut junctions_hq = JunctionStore::default();

        
        // Ensure cache is populated
        if read_provider.get_or_load_compact_cache(src).is_none() {
             populate_bam_cache(src, &progress, read_provider, &config.data_selection.genome_path)?;
        }

        // Tier 2: Retrieve and compute via compact reads (guaranteed to be there or fail above)
        let (compact_reads, segments, segment_tags, _names, ref_names, _ref_lengths, status) = read_provider.get_or_load_compact_cache(src)
            .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, format!("Failed to retrieve cache for {}", src)))?;
        
        if status != crate::reads_manager::CacheStatus::Complete {
             return Err(io::Error::new(io::ErrorKind::InvalidData, format!("Cache for {} is incomplete", src)));
        }

        println!("Rust: Calculating coverage from compact reads for {} ({} reads)", src, compact_reads.len());
        
        // We use the reference names from the cache, so we don't necessarily need the BAM header anymore.
        let mut coverage = Coverage::new(ref_names.len());
        
        for read in &compact_reads {
            coverage.add_compact_read(
                read, 
                &segments, 
                &segment_tags,
                &mut junctions, 
                &mut junctions_hq, 
                config.coverage.min_mapping_quality, 
                config.coverage.ambiguity.ambiguity_min_mapping_quality
            );
        }
        
        let mut sample_result = SampleResult {
            reference: "".to_string(),
            positions: Vec::new(),
            depths: Vec::new(),
            depths_hq: Vec::new(),
            depths_ambiguity: Vec::new(),
            junction_points: Vec::new(),
            junction_spans: Vec::new(),
            min_x: 0,
            max_x: 0,
            coverage_memory_mb: 0.0,
        };
        
        let local_covered_ranges = coverage.finalize(&ref_names, &junctions, &junctions_hq, &mut sample_result);
        sample_result.coverage_memory_mb = coverage.size_bytes() as f64 / (1024.0 * 1024.0);

        if config.data_selection.filter_outliers {
            sample_result.filter_outliers();
        } else {
            sample_result.detect_outliers();
        }

        println!("Rust: Finished processing sample {} in {:?}", sample_name, start_time.elapsed());
        
        let elapsed = start_time.elapsed();
        read_provider.benchmark.record_stage(&format!("Cohort: {}", sample_name), elapsed);
        read_provider.benchmark.observe_cohort_peak(&sample_name);

        Ok((sample_name, sample_result, local_covered_ranges))
    }).collect();

    let t_merging = Instant::now();
    let mut samples = HashMap::new();
    let mut global_min_x = usize::MAX;
    let mut global_max_x = 0;
    let mut cache_related_mb = 0.0;

    for res in results {
        let (name, s_res, l_ranges) = res?;
        cache_related_mb += s_res.coverage_memory_mb;
        
        if s_res.max_x > s_res.min_x {
            global_min_x = global_min_x.min(s_res.min_x);
            global_max_x = global_max_x.max(s_res.max_x);
        }

        samples.insert(name, s_res);
        for (ref_name, range) in l_ranges {
            let entry = all_covered_ranges.entry(ref_name).or_insert(range);
            entry.0 = entry.0.min(range.0);
            entry.1 = entry.1.max(range.1);
        }
    }
    
    read_provider.benchmark.record_stage("Result Merging", t_merging.elapsed());
    read_provider.benchmark.observe_peak();

    if global_min_x == usize::MAX {

        global_min_x = 0;
    }

    if !config.data_selection.filter_annotations {
        if let (Some(s), Some(e)) = (config.data_selection.analysis_start, config.data_selection.analysis_end) {
            println!("Rust: Using explicit metadata range for ROI focusing: {}..{}", s, e);
            global_min_x = s;
            global_max_x = e;
        }
    }

    let t_assembly = Instant::now();
    let assembly_report = config.data_selection.assembly_report_path.as_ref().and_then(|p| {
        let r = load_assembly_report(p).ok();
        if let Some(ref report) = r {
            println!("Rust: Loaded assembly report with {} mappings.", report.mapping.len());
        } else {
            println!("Rust: Failed to load assembly report from {}", p);
        }
        r
    });
    read_provider.benchmark.record_stage("Load Assembly Report", t_assembly.elapsed());
    read_provider.benchmark.observe_peak();
    
    progress.stage.store(2, Ordering::SeqCst); // 2: GTF Processing
    println!("Loading annotations from GTF...");
    
    let core_ranges = if !config.data_selection.filter_annotations {
        let mut ar = HashMap::new();
        for (ref_name, _) in &all_covered_ranges {
            ar.insert(ref_name.clone(), (global_min_x, global_max_x));
        }
        ar
    } else {
        all_covered_ranges.clone()
    };

    let padding = 500_000;
    let mut fetch_ranges = HashMap::new();
    for (ref_name, &(start, end)) in &core_ranges {
        let fstart = start.saturating_sub(padding);
        let fend = end.saturating_add(padding);
        fetch_ranges.insert(ref_name.clone(), (fstart, fend));
    }

    let mut annotations = Vec::new();
    let paths_to_process = if !config.data_selection.selected_gtfs.is_empty() {
        config.data_selection.selected_gtfs.clone()
    } else {
        let mut p = Vec::new();
        if let Some(first) = config.data_selection.gtf_paths.first() {
            p.push(first.clone());
        } else if let Some(path) = &config.data_selection.gtf_path {
            p.push(path.clone());
        }
        p
    };

    let t_gtf_total = Instant::now();
    for path in paths_to_process {
        let t_path = Instant::now();
        println!("Rust Debug: Calling get_annotations for path: {}", path);
        let filename = std::path::Path::new(&path)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or(&path);
        let offset = config.data_selection.gtf_offsets.get(filename).cloned().unwrap_or(0);
        
        let ann = get_annotations(&path, &core_ranges, &fetch_ranges, assembly_report.as_ref(), &progress, offset)?;
        println!("Rust Debug: Path {} returned {} annotations", path, ann.len());
        annotations.extend(ann);
        read_provider.benchmark.record_stage(&format!("GTF Load: {}", filename), t_path.elapsed());
        read_provider.benchmark.observe_peak();

        if annotations.len() >= 25_000 { break; }
    }
    read_provider.benchmark.record_stage("Total GTF Fetching", t_gtf_total.elapsed());

    
    println!("Found {} total annotations", annotations.len());

    progress.stage.store(3, Ordering::SeqCst); // 3: Finalizing JSON
    println!("Analysis complete. Returning to Python.");
    
    let t_final = Instant::now();

    let cache_core_mb = read_provider.get_core_cache_usage_mb();
    let cache_annotation_mb = read_provider.get_annotation_cache_usage_mb();

    let res = AnalysisResult { 
        samples, 
        annotations, 
        min_x: global_min_x, 
        max_x: global_max_x,
        cache_core_mb,
        cache_related_mb,
        cache_annotation_mb
    };
    read_provider.benchmark.record_stage("Pipeline Cleanup", t_final.elapsed());
    read_provider.benchmark.observe_peak();

    Ok(res)
}

