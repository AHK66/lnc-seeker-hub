// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

use lnc_seeker_lib::pipeline_config::AnyResult;
use lnc_seeker_lib::regions::RegionSet;
use noodles::core::Region;
use noodles::sam;
use sam::alignment::Record as _;
use crate::header_map::build_index_to_name;
use sam::alignment::RecordBuf;
use lnc_seeker_lib::reads_manager::IncrementalCacheBuilder;
use lnc_seeker_lib::genome::GenomeProvider;
use lnc_seeker_lib::models::ReadInfo;

/// Remap `rb`'s `reference_sequence_id` and `mate_reference_sequence_id`
/// from the source BAM's numeric indices (via `src_index_to_name`) into
/// the output header index space (`out_name_to_index`). Returns `true` if
/// the record is valid and should be kept.
fn remap_record_reference_ids(
    rb: &mut RecordBuf,
    bam_path: &str,
    src_index_to_name: &HashMap<usize, String>,
    out_name_to_index: &HashMap<String, usize>,
) -> bool {
    if let Some(src_rid) = rb.reference_sequence_id() {
        match src_index_to_name.get(&src_rid) {
            Some(src_name) => {
                if let Some(&out_rid) = out_name_to_index.get(src_name) {
                    *rb.reference_sequence_id_mut() = Some(out_rid);
                } else {
                    let q = rb
                        .name()
                        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                        .unwrap_or_else(|| "<unnamed>".to_string());
                    eprintln!(
                        "Skipping record {}: reference '{}' not present in output header (from {})",
                        q, src_name, bam_path
                    );
                    return false;
                }
            }
            None => {
                let q = rb
                    .name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                    .unwrap_or_else(|| "<unnamed>".to_string());
                eprintln!(
                    "Skipping record {}: invalid source ref id {} (from {})",
                    q, src_rid, bam_path
                );
                return false;
            }
        }
    }

    if let Some(src_mrid) = rb.mate_reference_sequence_id() {
        if let Some(mname) = src_index_to_name.get(&src_mrid) {
            if let Some(&out_mrid) = out_name_to_index.get(mname) {
                *rb.mate_reference_sequence_id_mut() = Some(out_mrid);
            } else {
                let q = rb
                    .name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                    .unwrap_or_else(|| "<unnamed>".to_string());
                eprintln!(
                    "Clearing mate ref for {}: mate reference '{}' not in output header (from {})",
                    q, mname, bam_path
                );
                *rb.mate_reference_sequence_id_mut() = None;
            }
        } else {
            *rb.mate_reference_sequence_id_mut() = None;
        }
    }

    true
}

/// Extract, remap and deduplicate `RecordBuf` alignments from the list of
/// `bam_paths` according to `remapped_introns`. Uses `out_name_to_index` to
/// remap reference IDs into the output header namespace. Returns collected
/// records grouped by region (inner vector index matches region index).
pub fn collect_records_from_bams(
    bam_paths: &[String],
    region_set: &RegionSet,
    annotate_endpoint_tag: bool,
    skip_secondary: bool,
    out_name_to_index: &HashMap<String, usize>,
    assembly_report: &lnc_seeker_lib::AssemblyReport,
    cohort_name: &str,
) -> AnyResult<Vec<Vec<RecordBuf>>> {
    let num_regions = region_set.regions().len();
    let mut collected: Vec<Vec<RecordBuf>> = vec![Vec::new(); num_regions];
    let mut collected_seen: Vec<HashSet<(String, Option<usize>, u16)>> = vec![HashSet::new(); num_regions];
    let mut skipped_duplicates: usize = 0;
    let mut skipped_examples: Vec<String> = Vec::new();

    let provider = lnc_seeker_lib::reads_manager::get_read_provider();

    for bam_path in bam_paths.iter() {
        println!("Processing BAM: {}", bam_path);
        provider.benchmark.observe_cohort_peak(cohort_name);

        let mut src_ireader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .map_err(|e| format!("{}: {}", bam_path, e))?;
        let src_header = src_ireader
            .read_header()
            .map_err(|e| format!("{}: {}", bam_path, e))?;
        let src_index_to_name = build_index_to_name(&src_header);
        // Build a set of reference names present in this BAM header so we can
        // avoid requesting regions that do not exist in the file (which
        // causes an error from the indexed reader).
        let src_names: HashSet<String> = src_header
            .reference_sequences()
            .iter()
            .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).into_owned())
            .collect();

        for (i, region) in region_set.regions().iter().enumerate() {
            let chrom = region.name().to_string();
            
            // Resolve which name is actually present in this BAM file
            let target_chrom = if src_names.contains(&chrom) {
                Some(chrom.clone())
            } else {
                assembly_report.reverse_mapping.get(&chrom)
                    .and_then(|cands| cands.iter().find(|&c| src_names.contains(c)).cloned())
            };

            let Some(actual_chrom) = target_chrom else {
                continue; // This BAM really doesn't have this chromosome
            };

            let query_region = if actual_chrom != chrom {
                Region::new(actual_chrom, region.interval())
            } else {
                region.clone()
            };

            provider.benchmark.observe_cohort_peak(cohort_name);

            let min_start = query_region.interval().start().map(usize::from).unwrap_or(1);
            let max_end = query_region.interval().end().map(usize::from).unwrap_or(usize::MAX);

            for result in src_ireader
                .query(&src_header, &query_region)
                .map_err(|e| format!("{}: {}: {}", bam_path, query_region, e))?
                .records()
            {
                let record = result.map_err(|e| format!("{}: {}", bam_path, e))?;

                if skip_secondary && record.flags().is_secondary() {
                    continue;
                }

                let aln_start_opt = record
                    .alignment_start()
                    .transpose()
                    .map_err(|e| format!("{}: {}", bam_path, e))?;
                if aln_start_opt.is_none() {
                    continue;
                }

                let s = usize::from(aln_start_opt.unwrap());
                let alignment_span_opt = record
                    .alignment_span()
                    .transpose()
                    .map_err(|e| format!("{}: {}", bam_path, e))?;
                let e = if let Some(span) = alignment_span_opt { if span > 0 { s + span - 1 } else { s } } else { s };

                // Fix: Correct overlap check to include reads that span the entire region
                let overlaps = s <= max_end && e >= min_start;

                if !overlaps {
                    continue;
                }

                let mut rb = RecordBuf::try_from_alignment_record(&src_header, &record)
                    .map_err(|e| format!("{}: {}", bam_path, e))?;

                // Remap reference ids into output header namespace.
                if !remap_record_reference_ids(&mut rb, bam_path, &src_index_to_name, out_name_to_index) {
                    continue;
                }

                // Optionally attach endpoint tag
                if annotate_endpoint_tag {
                    let start_in = s >= min_start && s <= max_end;
                    let end_in = e >= min_start && e <= max_end;
                    let indicator = match (start_in, end_in) {
                        (true, true) => "LR",
                        (true, false) => "L",
                        (false, true) => "R",
                        _ => "",
                    };
                    if !indicator.is_empty() {
                        rb.data_mut().insert(
                            noodles::sam::alignment::record::data::field::Tag::new(b'E', b'P'),
                            noodles::sam::alignment::record_buf::data::field::Value::from(indicator),
                        );
                    }
                }

                let base_name = rb
                    .name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                    .unwrap_or_else(|| "unnamed".to_string());
                let start_opt = rb.alignment_start().map(|p| usize::from(p));
                let flags = u16::from(rb.flags());
                let key = (base_name.clone(), start_opt, flags);
                if collected_seen[i].insert(key) {
                    collected[i].push(rb);
                } else {
                    skipped_duplicates += 1;
                    if skipped_examples.len() < 10 {
                        let start_str = start_opt
                            .map(|s| s.to_string())
                            .unwrap_or_else(|| "unknown".to_string());
                        skipped_examples.push(format!("{} start={} flags={:04x} (from {})", base_name, start_str, flags, bam_path));
                    }
                }
            }
        }
    }

    if skipped_duplicates > 0 {
        eprintln!(
            "Skipped {} duplicate collected records (showing up to {} examples):",
            skipped_duplicates,
            skipped_examples.len()
        );
        for ex in skipped_examples.iter() {
            eprintln!("  {}", ex);
        }
    }

    Ok(collected)
}

/// Streamed version of record collection that avoids large memory spikes by
/// converting to compact storage incrementally.
pub fn collect_records_streamed(
    bam_paths: &[String],
    region_set: &RegionSet,
    skip_secondary: bool,
    out_header: &sam::Header,
    assembly_report: &lnc_seeker_lib::AssemblyReport,
    genome_provider: &mut GenomeProvider,
    cohort_name: &str,
) -> AnyResult<Vec<lnc_seeker_lib::reads_manager::BamCache>> {
    let t_start = std::time::Instant::now();
    let num_regions = region_set.regions().len();
    let provider = lnc_seeker_lib::reads_manager::get_read_provider();
    
    let mut ref_names = Vec::new();
    let mut ref_lengths = Vec::new();
    for (name, seq) in out_header.reference_sequences().iter() {
        let name_bytes: &[u8] = name.as_ref();
        ref_names.push(String::from_utf8_lossy(name_bytes).to_string());
        ref_lengths.push(seq.length().get());
    }

    let mut builders: Vec<IncrementalCacheBuilder> = Vec::with_capacity(num_regions);
    for j in 0..num_regions {
        let region_label = region_set.regions()[j].to_string();
        builders.push(IncrementalCacheBuilder::new(&region_label, ref_names.clone(), ref_lengths.clone()));
    }

    let mut collected_seen: Vec<HashSet<u64>> = vec![HashSet::new(); num_regions];
    let mut skipped_duplicates: usize = 0;

    for bam_path in bam_paths.iter() {
        println!("Processing BAM (Streamed): {}", bam_path);
        provider.benchmark.observe_cohort_peak(cohort_name);

        let mut src_ireader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .map_err(|e| format!("{}: {}", bam_path, e))?;
        let src_header = src_ireader
            .read_header()
            .map_err(|e| format!("{}: {}", bam_path, e))?;
        
        let src_names: HashSet<String> = src_header
            .reference_sequences()
            .iter()
            .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).into_owned())
            .collect();

        for (i, region) in region_set.regions().iter().enumerate() {
            let chrom = region.name().to_string();
            
            let target_chrom = if src_names.contains(&chrom) {
                Some(chrom.clone())
            } else {
                assembly_report.reverse_mapping.get(&chrom)
                    .and_then(|cands| cands.iter().find(|&c| src_names.contains(c)).cloned())
            };

            let Some(actual_chrom) = target_chrom else {
                continue;
            };

            let query_region = if actual_chrom != chrom {
                Region::new(actual_chrom.clone(), region.interval())
            } else {
                region.clone()
            };

            provider.benchmark.observe_cohort_peak(cohort_name);

            let min_start = query_region.interval().start().map(usize::from).unwrap_or(1);
            let max_end = query_region.interval().end().map(usize::from).unwrap_or(usize::MAX);

            // Pre-fetch reference sequence for mismatch calculation
            let ref_seq_full = genome_provider.get_sequence(&actual_chrom).ok().flatten();

            for result in src_ireader
                .query(&src_header, &query_region)
                .map_err(|e| format!("{}: {}: {}", bam_path, query_region, e))?
                .records()
            {
                let record = result.map_err(|e| format!("{}: {}", bam_path, e))?;

                if skip_secondary && record.flags().is_secondary() {
                    continue;
                }

                let aln_start_opt = record
                    .alignment_start()
                    .transpose()
                    .map_err(|e| format!("{}: {}", bam_path, e))?;
                if aln_start_opt.is_none() {
                    continue;
                }

                let s = usize::from(aln_start_opt.unwrap());
                let alignment_span_opt = record
                    .alignment_span()
                    .transpose()
                    .map_err(|e| format!("{}: {}", bam_path, e))?;
                let e = if let Some(span) = alignment_span_opt { if span > 0 { s + span - 1 } else { s } } else { s };

                if s > max_end || e < min_start {
                    continue;
                }

                // Deduplicate using hash
                let mut hasher = std::collections::hash_map::DefaultHasher::new();
                let name_bytes: &[u8] = record.name().map(|n| n.as_ref()).unwrap_or(b"unnamed");
                name_bytes.hash(&mut hasher);
                s.hash(&mut hasher);
                u16::from(record.flags()).hash(&mut hasher);
                let h = hasher.finish();

                if !collected_seen[i].insert(h) {
                    skipped_duplicates += 1;
                    continue;
                }

                // Conversion to ReadInfo
                let ref_seq_slice = ref_seq_full.as_ref().map(|s| s.as_ref());
                if let Some(mut info) = ReadInfo::from_record(&record, &src_header, ref_seq_slice)? {
                    // Normalize reference name to the one expected in out_header
                    info.reference = chrom.clone();
                    builders[i].add_read(info);
                }
            }
        }
    }

    if skipped_duplicates > 0 {
        println!("Deduplication (Streamed): Skipped {} duplicate records.", skipped_duplicates);
    }

    let results = builders.into_iter()
        .map(|b| {
            let bam_path = b.bam_path.clone();
            let mapping_mb = b.get_mapping_memory_estimate() as f64 / 1024.0 / 1024.0;
            let mut cache = b.finalize(lnc_seeker_lib::reads_manager::CacheStatus::Complete);
            
            // Record stats before compression
            let uncompressed_header = cache.get_uncompressed_header_bytes();
            let payload = cache.get_payload_bytes();
            let read_count = cache.reads.len();
            let segment_count = cache.segments.len();
            let data_mb = cache.size_bytes as f64 / 1024.0 / 1024.0;
            
            // Compress headers immediately to lower memory footprint before returning
            cache.compress_headers();
            
            let provider = lnc_seeker_lib::reads_manager::get_read_provider();
            provider.benchmark.record_bam_stats(lnc_seeker_lib::benchmarking::BamStats {
                context: bam_path,
                read_count,
                segment_count,
                data_mb,
                mapping_mb,
                os_mb: lnc_seeker_lib::benchmarking::get_os_memory() as f64 / 1024.0 / 1024.0,
                header_uncompressed_bytes: uncompressed_header,
                header_compressed_bytes: cache.get_compressed_header_bytes(),
                payload_uncompressed_bytes: payload,
                payload_compressed_bytes: payload,
            });

            cache
        })
        .collect();

    let provider = lnc_seeker_lib::reads_manager::get_read_provider();
    provider.benchmark.record_stage("Streamed Collection (All Regions)", t_start.elapsed());

    Ok(results)
}

// Summarize duplicate skipping when collection finishes.
// Note: keep this function lightweight; the caller prints the returned
// records. We print a short summary on stderr so users see the reduction
// in noise without losing visibility into examples.
// (Placed after the function to avoid breaking return flow in older callers.)
