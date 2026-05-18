// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::collections::{HashMap, HashSet};
use noodles::core::Region;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// Centralized handling of NCBI assembly reports for contig name translation.
#[derive(Debug, Clone)]
pub struct AssemblyReport {
    pub mapping: HashMap<String, String>,        // Any identifier -> UCSC-style name
    pub reverse_mapping: HashMap<String, Vec<String>>, // UCSC-style name -> [all identifiers]
    pub mapped_names: HashSet<String>,           // Set of known UCSC-style names
}

impl AssemblyReport {
    /// Parse an NCBI assembly report file.
    /// Maps identifiers (Sequence-Name, GenBank-Accn, RefSeq-Accn) 
    /// to UCSC-style display names (column 10).
    pub fn parse_from_file<P: AsRef<std::path::Path>>(path: P) -> io::Result<Self> {
        let f = File::open(path)?;
        let mut mapping: HashMap<String, String> = HashMap::new();
        let mut reverse_mapping: HashMap<String, Vec<String>> = HashMap::new();
        let mut mapped_names: HashSet<String> = HashSet::new();

        for line in BufReader::new(f).lines() {
            let l = line?;
            if l.trim().is_empty() || l.starts_with('#') {
                continue;
            }

            let cols: Vec<&str> = l.split('\t').collect();
            // Need at least 10 columns for UCSC-style-name at index 9
            if cols.len() < 10 {
                continue;
            }
            
            let ucsc = cols[9].trim().to_string();
            if ucsc.is_empty() || ucsc == "-" {
                continue;
            }

            mapped_names.insert(ucsc.clone());
            
            // Collect all potential names for this chromosome
            let candidates = vec![
                cols[0].trim().to_string(), // Sequence-Name
                cols[4].trim().to_string(), // GenBank-Accn
                cols[6].trim().to_string(), // RefSeq-Accn
            ];

            for cand in candidates {
                if !cand.is_empty() && cand != "-" {
                    mapping.insert(cand.clone(), ucsc.clone());
                    reverse_mapping.entry(ucsc.clone()).or_default().push(cand);
                }
            }
        }
        Ok(Self { mapping, reverse_mapping, mapped_names })
    }

    /// Translates a region's name if a mapping exists.
    pub fn translate_region(&self, region: &Region) -> Region {
        let name = String::from_utf8_lossy(region.name().as_ref()).into_owned();
        if self.mapped_names.contains(&name) {
            return region.clone();
        }
        if let Some(mapped) = self.mapping.get(&name) {
            return Region::new(mapped.clone(), region.interval());
        }
        region.clone()
    }
}
