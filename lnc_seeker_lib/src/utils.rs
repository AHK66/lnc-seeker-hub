// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
/// Shared attributes for genomic features (GTF/GFF).
#[derive(Debug, Clone, Default)]
pub struct FeatureAttributes {
    pub gene_id: String,
    pub gene_name: Option<String>,
    pub transcript_id: String,
    pub exon_number: Option<String>,
}

/// Parses the attributes string from a GTF/GFF line.
/// Supports both GTF (space-delimited) and GFF3 (equals-delimited) styles.
pub fn parse_attributes(attribute_str: &str) -> FeatureAttributes {
    let mut attr = FeatureAttributes::default();
    attr.gene_id = "unknown".to_string();
    attr.transcript_id = "unknown".to_string();

    for part in attribute_str.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }

        if part.contains('=') {
            // GFF3 style
            let kv: Vec<&str> = part.splitn(2, '=').collect();
            if kv.len() == 2 {
                let k = kv[0].trim().to_lowercase();
                let v = kv[1]
                    .trim()
                    .trim_matches('"')
                    .trim_matches('\'')
                    .to_string();
                match k.as_str() {
                    "id" => {
                        if attr.transcript_id == "unknown" {
                            attr.transcript_id = v;
                        }
                    }
                    "parent" => {
                        if attr.gene_id == "unknown" {
                            attr.gene_id = v;
                        }
                    }
                    "name" | "gene_name" | "gene" => {
                        if attr.gene_name.is_none() {
                            attr.gene_name = Some(v);
                        }
                    }
                    "gene_id" => {
                        if attr.gene_id == "unknown" {
                            attr.gene_id = v;
                        }
                    }
                    "transcript_id" => {
                        if attr.transcript_id == "unknown" {
                            attr.transcript_id = v;
                        }
                    }
                    _ => {}
                }
            }
        } else {
            // GTF style: key "value"
            let mut items = part.split_whitespace();
            if let Some(key) = items.next() {
                let val = items
                    .collect::<Vec<_>>()
                    .join(" ")
                    .trim_matches('"')
                    .trim_matches('\'')
                    .to_string();
                match key {
                    "gene_id" => attr.gene_id = val,
                    "gene_name" | "gene" => attr.gene_name = Some(val),
                    "transcript_id" => attr.transcript_id = val,
                    "exon_number" => attr.exon_number = Some(val),
                    _ => {}
                }
            }
        }
    }
    attr
}

/// Checks if a record at [start, end] overlaps with a range [q_start, q_end].
pub fn intersects(start: usize, end: usize, q_start: usize, q_end: usize) -> bool {
    start < q_end && end > q_start
}

pub fn normalize_path(path: &str) -> String {
    let p = path.replace('\\', "/");
    #[cfg(windows)]
    {
        p.to_lowercase()
    }
    #[cfg(not(windows))]
    {
        p
    }
}
