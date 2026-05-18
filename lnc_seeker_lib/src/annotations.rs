// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufReader, BufRead};
use std::sync::atomic::Ordering;

use noodles::csi;
use flate2::read::MultiGzDecoder;
use noodles::core::Region;

use crate::models::{Annotation, PackedAnnotation, LbaFile};
use crate::progress::ProgressData;
use crate::assembly::AssemblyReport;
use crate::utils::{parse_attributes, intersects};

pub fn get_annotations(
    path: &str,
    core_ranges: &HashMap<String, (usize, usize)>,
    fetch_ranges: &HashMap<String, (usize, usize)>,
    report_opt: Option<&AssemblyReport>,
    progress: &ProgressData,
    offset: i32,
) -> io::Result<Vec<Annotation>> {
    // Check if optimized binary version exists
    let lba_path = format!("{}.lba", path);
    if std::path::Path::new(&lba_path).exists() {
        println!("Rust: Found optimized LBA file: {}", lba_path);
        match get_annotations_lba(&lba_path, core_ranges, fetch_ranges, report_opt, offset) {
            Ok(res) => {
                println!("Rust: Successfully loaded {} features from LBA", res.len());
                return Ok(res);
            }
            Err(e) => {
                println!("Warning: Failed to load LBA file (falling back to GTF): {}", e);
            }
        }
    }

    println!("Rust: Attempting indexed GTF access for: {}", path);
    let mut results: Vec<Annotation> = Vec::new();

    let index_path_candidates = vec![format!("{}.tbi", path), format!("{}.csi", path)];
    let mut index_opt = None;
    for p in index_path_candidates {
        if std::path::Path::new(&p).exists() {
            if let Ok(idx) = csi::fs::read(&p) {
                index_opt = Some(idx);
                break;
            }
        }
    }

    if let Some(index) = index_opt {
        println!("Rust: Opened CSI/Tabix index for: {}", path);
        let mut ireader = csi::io::IndexedReader::new(File::open(path)?, index);

        // --- PASS 1: Discovery ---
        let mut discovered_transcripts = HashSet::new();
        for (norm_ref, &(rstart, rend)) in core_ranges.iter() {
            if rstart >= rend { continue; }
            let candidate_seqs = get_candidate_seq_names(norm_ref, report_opt);
            for seq_name in candidate_seqs.iter() {
                let region_str = format!("{}:{}-{}", seq_name, rstart + 1, rend);
                let region: Region = match region_str.parse() {
                    Ok(r) => r,
                    Err(_) => continue,
                };
                let query_iter = match ireader.query(&region) {
                    Ok(it) => it,
                    Err(_) => continue,
                };
                for rec_result in query_iter {
                    let rec = rec_result?;
                    let line = rec.as_ref();
                    if line.starts_with('#') || line.is_empty() { continue; }
                    
                    // Optimization: Use nth(8) to get attribute field without collecting into Vec
                    if let Some(attr_part) = line.split('\t').nth(8) {
                        let attr = parse_attributes(attr_part);
                        if !attr.transcript_id.is_empty() {
                            discovered_transcripts.insert(attr.transcript_id);
                        }
                    }
                }
            }
        }
        println!("Rust: Pass 1 discovered {} transcripts in ROI", discovered_transcripts.len());

        // --- PASS 2: Collection ---
        let mut total_found = 0usize;
        for (norm_ref, &(fstart, fend)) in fetch_ranges.iter() {
            if fstart >= fend { continue; }
            let candidate_seqs = get_candidate_seq_names(norm_ref, report_opt);
            for seq_name in candidate_seqs.iter() {
                let region_str = format!("{}:{}-{}", seq_name, fstart + 1, fend);
                let region: Region = match region_str.parse() {
                    Ok(r) => r,
                    Err(_) => continue,
                };
                let query_iter = match ireader.query(&region) {
                    Ok(it) => it,
                    Err(_) => continue,
                };
                for rec_result in query_iter {
                    let rec = rec_result?;
                    let line = rec.as_ref();
                    if line.starts_with('#') || line.is_empty() { continue; }
                    
                    // Optimization: Access fields via iterator to avoid allocations
                    let mut fields = line.split('\t');
                    let chrom = fields.next().unwrap_or("");
                    let _source = fields.next();
                    let feature_type = fields.next().unwrap_or("");
                    let start_str = fields.next().unwrap_or("");
                    let end_str = fields.next().unwrap_or("");
                    let _score = fields.next();
                    let strand = fields.next().unwrap_or("");
                    let _frame = fields.next();
                    let attr_str = fields.next().unwrap_or("");

                    if !is_relevant_feature(feature_type) { continue; }
                    let attr = parse_attributes(attr_str);
                    if !discovered_transcripts.contains(&attr.transcript_id) { continue; }

                    let raw_start = start_str.parse::<i64>().unwrap_or(1);
                    let start = (raw_start + offset as i64).saturating_sub(1).max(0) as usize;
                    let end = end_str.parse::<usize>().unwrap_or(0);

                    let normalized = if let Some(rep) = report_opt {
                        rep.mapping.get(chrom).cloned().unwrap_or_else(|| chrom.to_string())
                    } else {
                        chrom.to_string()
                    };

                    results.push(Annotation {
                        reference: normalized.clone(),
                        start,
                        end,
                        feature: feature_type.to_string(),
                        gene_id: attr.gene_id,
                        gene_name: attr.gene_name,
                        transcript_id: attr.transcript_id,
                        exon_number: attr.exon_number,
                        strand: strand.to_string(),
                    });
                    total_found += 1;
                    if total_found >= 25_000 {
                        println!("Rust: Reached limit of 25,000 annotations.");
                        return Ok(results);
                    }
                }
            }
        }
        return Ok(results);
    }

    // Streaming fallback (Single-pass over core_ranges)
    println!("Rust: Indexed reader unavailable — falling back to streaming parser for: {}", path);
    let file = File::open(path)?;
    let decoder = MultiGzDecoder::new(file);
    let reader = BufReader::new(decoder);
    let mut total_found = 0usize;
    for (count, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        if line.starts_with('#') || line.is_empty() { continue; }
        
        let mut fields = line.split('\t');
        let chrom = fields.next().unwrap_or("");
        let _source = fields.next();
        let feature_type = fields.next().unwrap_or("");
        let start_str = fields.next().unwrap_or("");
        let end_str = fields.next().unwrap_or("");
        let _score = fields.next();
        let strand = fields.next().unwrap_or("");
        let _frame = fields.next();
        let attr_str = fields.next().unwrap_or("");

        if !is_relevant_feature(feature_type) { continue; }
        
        let normalized = if let Some(rep) = report_opt {
            rep.mapping.get(chrom).cloned().unwrap_or_else(|| chrom.to_string())
        } else {
            chrom.to_string()
        };
        
        if let Some(&(rstart, rend)) = core_ranges.get(&normalized) {
            let raw_start = start_str.parse::<i64>().unwrap_or(1);
            let start = (raw_start + offset as i64).saturating_sub(1).max(0) as usize;
            let end = end_str.parse::<usize>().unwrap_or(0);
            if intersects(start, end, rstart, rend) {
                let attr = parse_attributes(attr_str);
                results.push(Annotation {
                    reference: normalized.clone(),
                    start,
                    end,
                    feature: feature_type.to_string(),
                    gene_id: attr.gene_id,
                    gene_name: attr.gene_name,
                    transcript_id: attr.transcript_id,
                    exon_number: attr.exon_number,
                    strand: strand.to_string(),
                });
                total_found += 1;
                if total_found >= 25_000 { break; }
            }
        }
        if count % 100_000 == 0 { progress.current.store(count as u32, Ordering::Relaxed); }
    }
    Ok(results)
}

fn get_candidate_seq_names(norm_ref: &str, report_opt: Option<&AssemblyReport>) -> Vec<String> {
    let mut candidate_seqs: Vec<String> = Vec::new();
    candidate_seqs.push(norm_ref.to_string());
    if let Some(report) = report_opt {
        if let Some(ucsc) = report.mapping.get(norm_ref) {
            candidate_seqs.push(ucsc.clone());
            if let Some(others) = report.reverse_mapping.get(ucsc) {
                for o in others { candidate_seqs.push(o.clone()); }
            }
        } else if report.mapped_names.contains(norm_ref) {
            if let Some(others) = report.reverse_mapping.get(norm_ref) {
                for o in others { candidate_seqs.push(o.clone()); }
            }
        }
    }
    if norm_ref.starts_with("chr") {
        candidate_seqs.push(norm_ref.replace("chr", ""));
    } else if norm_ref.len() <= 2 || (norm_ref == "M" || norm_ref == "MT") {
        candidate_seqs.push(format!("chr{}", norm_ref));
    }
    let mut seen_ids = HashSet::new();
    candidate_seqs.retain(|s| seen_ids.insert(s.clone()));
    candidate_seqs
}

fn is_relevant_feature(f: &str) -> bool {
    matches!(f, "exon" | "CDS" | "five_prime_utr" | "three_prime_utr" | "stop_codon" | "start_codon" | "ncRNA" | "transcript" | "non_coding_exon" | "lnc_RNA")
}

/// Converts a GTF file to the optimized binary LBA (LNC-Seeker Binary Annotation) format.
pub fn optimize_gtf_to_lba(gtf_path: &str, lba_path: &str) -> io::Result<()> {
    println!("LBA: Optimizing {} -> {}...", gtf_path, lba_path);
    let file = File::open(gtf_path)?;
    let decoder = MultiGzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let mut chroms = Vec::new();
    let mut chrom_to_idx = HashMap::new();
    let mut string_pool = Vec::new();
    let mut string_to_idx = HashMap::new();
    let mut records = Vec::new();

    let get_string_idx = |string_to_idx: &mut HashMap<String, u32>, string_pool: &mut Vec<String>, s: String| -> u32 {
        if let Some(&idx) = string_to_idx.get(&s) {
            idx
        } else {
            let idx = string_pool.len() as u32;
            string_to_idx.insert(s.clone(), idx);
            string_pool.push(s);
            idx
        }
    };

    let get_chrom_idx = |chrom_to_idx: &mut HashMap<String, u16>, chroms: &mut Vec<String>, s: &str| -> u16 {
        if let Some(&idx) = chrom_to_idx.get(s) {
            idx
        } else {
            let idx = chroms.len() as u16;
            chrom_to_idx.insert(s.to_string(), idx);
            chroms.push(s.to_string());
            idx
        }
    };

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with('#') || line.is_empty() { continue; }
        
        let mut fields = line.split('\t');
        let chrom = fields.next().unwrap_or("");
        let _source = fields.next();
        let feature_type = fields.next().unwrap_or("");
        let start_str = fields.next().unwrap_or("");
        let end_str = fields.next().unwrap_or("");
        let _score = fields.next();
        let strand = fields.next().unwrap_or("");
        let _frame = fields.next();
        let attr_str = fields.next().unwrap_or("");

        if !is_relevant_feature(feature_type) { continue; }

        let attr = parse_attributes(attr_str);
        
        // Compact mapping for common features (case-insensitive)
        let feature_idx = match feature_type.to_lowercase().as_str() {
            "exon" => 1,
            "transcript" => 2,
            "cds" => 3,
            "gene" => 4,
            "five_prime_utr" | "three_prime_utr" | "utr" => 5,
            _ => 0,
        };

        let strand_idx = match strand {
            "+" => 0,
            "-" => 1,
            _ => 2,
        };

        records.push(PackedAnnotation {
            chrom_idx: get_chrom_idx(&mut chrom_to_idx, &mut chroms, chrom),
            start: start_str.parse::<u32>().unwrap_or(0).saturating_sub(1),
            end: end_str.parse::<u32>().unwrap_or(0),
            feature_idx,
            strand: strand_idx,
            gene_id_idx: get_string_idx(&mut string_to_idx, &mut string_pool, attr.gene_id),
            gene_name_idx: attr.gene_name.map(|s| get_string_idx(&mut string_to_idx, &mut string_pool, s)).unwrap_or(u32::MAX),
            transcript_id_idx: get_string_idx(&mut string_to_idx, &mut string_pool, attr.transcript_id),
            exon_number_idx: attr.exon_number.map(|s| get_string_idx(&mut string_to_idx, &mut string_pool, s)).unwrap_or(u32::MAX),
        });
    }

    // Sort records by chrom_idx and then by start position for fast binary search
    records.sort_by(|a, b| {
        a.chrom_idx.cmp(&b.chrom_idx)
            .then(a.start.cmp(&b.start))
    });

    let lba = LbaFile {
        magic: *b"LBA\x01",
        chroms,
        string_pool,
        records,
    };

    let mut out = File::create(lba_path)?;
    bincode::serialize_into(&mut out, &lba).map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    
    println!("LBA: Successfully converted {} records to binary format.", lba.records.len());
    Ok(())
}

fn get_annotations_lba(
    path: &str,
    core_ranges: &HashMap<String, (usize, usize)>,
    fetch_ranges: &HashMap<String, (usize, usize)>,
    report_opt: Option<&AssemblyReport>,
    offset: i32,
) -> io::Result<Vec<Annotation>> {
    let lba = crate::reads_manager::get_read_provider().get_or_load_lba(path)?;
    
    // --- PASS 1: Discovery (Transcripts in ROI) ---
    let mut discovered_transcripts = HashSet::new();
    for (norm_ref, &(rstart, rend)) in core_ranges.iter() {
        let candidates = get_candidate_seq_names(norm_ref, report_opt);
        for name in &candidates {
            if let Some(c_idx) = lba.chroms.iter().position(|x| x == name) {
                let c_idx = c_idx as u16;
                // Robust slice search: Binary search for chrom start, then linear scan
                let start_idx = lba.records.partition_point(|r| r.chrom_idx < c_idx);
                
                for i in start_idx..lba.records.len() {
                    let r = &lba.records[i];
                    if r.chrom_idx != c_idx { break; }
                    
                    // Simple overlap check (0-based start, 1-based end)
                    if r.start < rend as u32 && r.end >= rstart as u32 {
                        let t_id = lba.string_pool.get(r.transcript_id_idx as usize).cloned().unwrap_or_default();
                        if !t_id.is_empty() {
                            discovered_transcripts.insert(t_id);
                        }
                    }
                }
            }
        }
    }

    // --- PASS 2: Collection (All features of discovered transcripts in fetch_ranges) ---
    let mut results = Vec::new();
    for (norm_ref, &(fstart, fend)) in fetch_ranges.iter() {
        let candidates = get_candidate_seq_names(norm_ref, report_opt);
        for name in &candidates {
            if let Some(c_idx) = lba.chroms.iter().position(|x| x == name) {
                let c_idx = c_idx as u16;
                let start_idx = lba.records.partition_point(|r| r.chrom_idx < c_idx);
                
                for i in start_idx..lba.records.len() {
                    let r = &lba.records[i];
                    if r.chrom_idx != c_idx { break; }
                    
                    if r.start < fend as u32 && r.end >= fstart as u32 {
                        let t_id = lba.string_pool.get(r.transcript_id_idx as usize).cloned().unwrap_or_default();
                        if discovered_transcripts.contains(&t_id) {
                            let feature = match r.feature_idx {
                                1 => "exon",
                                2 => "transcript",
                                3 => "CDS",
                                4 => "gene",
                                5 => "UTR",
                                _ => "unknown",
                            };

                            let start = (r.start as i64 + offset as i64).max(0) as usize;
                            let end = r.end as usize;

                            results.push(Annotation {
                                reference: norm_ref.clone(),
                                start,
                                end,
                                feature: feature.to_string(),
                                gene_id: lba.string_pool.get(r.gene_id_idx as usize).cloned().unwrap_or_default(),
                                gene_name: if r.gene_name_idx == u32::MAX { None } else { lba.string_pool.get(r.gene_name_idx as usize).cloned() },
                                transcript_id: t_id,
                                exon_number: if r.exon_number_idx == u32::MAX { None } else { lba.string_pool.get(r.exon_number_idx as usize).cloned() },
                                strand: match r.strand {
                                    0 => "+",
                                    1 => "-",
                                    _ => ".",
                                }.to_string(),
                            });
                        }
                    }
                    if results.len() >= 25_000 { return Ok(results); }
                }
            }
        }
    }
    
    Ok(results)
}
