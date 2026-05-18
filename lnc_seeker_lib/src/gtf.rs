// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
//! Queries a bgzip-compressed GTF file and maps gene names to regions.

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;

use noodles::bgzf;
use noodles::core::{Position, Region};
use noodles::gff;
use noodles::gtf;
use crate::AssemblyReport;

/// Extract all GTF records whose `gene_name` or `gene_id` attribute matches any of the `names`.
pub fn extract_gene_names(
    src: &str,
    names: &[String],
) -> Result<Vec<gff::feature::RecordBuf>, Box<dyn Error>> {
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(bgzf::io::Reader::new)
        .map(gtf::io::Reader::new)?;

    let mut matches: Vec<gff::feature::RecordBuf> = Vec::new();

    let names_lower: HashSet<String> = names.iter().map(|n| n.to_ascii_lowercase()).collect();

    for result in reader.record_bufs() {
        let record = result?;

        let mut matched = false;

        if let Some(val) = record.attributes().get(b"gene_name") {
            if let Some(bs) = val.as_string() {
                if let Ok(s) = std::str::from_utf8(bs.as_ref()) {
                    if names_lower.contains(&s.to_ascii_lowercase()) {
                        matched = true;
                    }
                }
            }
        }

        if !matched {
            if let Some(val) = record.attributes().get(b"gene_id") {
                if let Some(bs) = val.as_string() {
                    if let Ok(s) = std::str::from_utf8(bs.as_ref()) {
                        if names_lower.contains(&s.to_ascii_lowercase()) {
                            matched = true;
                        }
                    }
                }
            }
        }

        if matched {
            matches.push(record);
        }
    }

    Ok(matches)
}

/// Returns a map of gene names to their spanned genomic regions using an `AssemblyReport`.
pub fn get_gene_regions(
    gtf_path: &str,
    assembly_report: &AssemblyReport,
    gene_names: &[String],
) -> Result<HashMap<String, Region>, Box<dyn Error>> {
    let all_matches = extract_gene_names(gtf_path, gene_names)?;

    let mut gene_to_data: HashMap<String, (String, usize, usize)> = HashMap::new();

    let names_lower: HashSet<String> = gene_names.iter().map(|n| n.to_lowercase()).collect();

    for record in all_matches {
        let mut name = None;
        if let Some(val) = record.attributes().get(b"gene_name") {
            if let Some(bs) = val.as_string() {
                if let Ok(s) = std::str::from_utf8(bs.as_ref()) {
                    if names_lower.contains(&s.to_lowercase()) {
                        name = Some(s.to_string());
                    }
                }
            }
        }
        if name.is_none() {
            if let Some(val) = record.attributes().get(b"gene_id") {
                if let Some(bs) = val.as_string() {
                    if let Ok(s) = std::str::from_utf8(bs.as_ref()) {
                        if names_lower.contains(&s.to_lowercase()) {
                            name = Some(s.to_string());
                        }
                    }
                }
            }
        }

        let Some(gene_name) = name else { continue; };
        
        let canonical_name = gene_names.iter()
            .find(|&n| n.to_lowercase() == gene_name.to_lowercase())
            .cloned()
            .unwrap_or(gene_name);

        let chrom = record.reference_sequence_name();
        let start = usize::from(record.start());
        let end = usize::from(record.end());

        let mapped_chrom = assembly_report.mapping.get(chrom).cloned().unwrap_or_else(|| chrom.to_string());

        let entry = gene_to_data.entry(canonical_name).or_insert((mapped_chrom.clone(), usize::MAX, 0));
        entry.1 = entry.1.min(start);
        entry.2 = entry.2.max(end);
    }

    let mut result = HashMap::new();
    for (name, (chrom, start, end)) in gene_to_data {
        let start_pos = Position::try_from(start)?;
        let end_pos = Position::try_from(end)?;
        result.insert(name, Region::new(chrom, start_pos..=end_pos));
    }

    Ok(result)
}
