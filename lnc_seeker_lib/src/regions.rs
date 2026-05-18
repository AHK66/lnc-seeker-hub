// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use noodles::core::{Region, Position};
use std::collections::HashMap;
use std::convert::TryFrom;
use crate::AssemblyReport;

/// A set of genomic regions for read extraction.
#[derive(Debug, Clone)]
pub struct RegionSet {
    pub(crate) regions: Vec<Region>,
}

impl RegionSet {
    /// Creates a new `RegionSet` from a vector of `Region`s, checking and translating
    /// chromosome names using the provided `AssemblyReport`.
    pub fn new(regions: Vec<Region>, assembly_report: &AssemblyReport) -> Self {
        let regions = regions
            .into_iter()
            .map(|r| assembly_report.translate_region(&r))
            .collect();
        Self { regions }
    }

    /// Returns the regions in the set.
    pub fn regions(&self) -> &[Region] {
        &self.regions
    }

    /// Returns the number of regions in the set.
    pub fn len(&self) -> usize {
        self.regions.len()
    }

    /// Filters the regions in the set to only include those at the specified indices.
    pub fn filter_by_indices(&mut self, indices: &[usize]) {
        self.regions = indices.iter()
            .filter_map(|&i| self.regions.get(i))
            .cloned()
            .collect();
    }

    /// Pretty-prints the regions in the set to stdout.
    pub fn pretty_print(&self) {
        println!("Extraction Regions:");
        if self.regions.is_empty() {
            println!("  (none)");
        } else {
            for region in &self.regions {
                println!("  - {}", region);
            }
        }
    }
}

/// Returns the overall region spanned by the given genes, with an optional offset.
/// All genes must be on the same chromosome.
pub fn get_overall_gene_region(
    gene_names: Vec<String>,
    gene_regions: &HashMap<String, Region>,
    offset: usize,
) -> Result<Region, String> {
    let mut chrom: Option<String> = None;
    let mut min_start = usize::MAX;
    let mut max_end = 0;

    for name in gene_names {
        let region = gene_regions.get(&name).ok_or_else(|| format!("gene not found: {}", name))?;
        let c = String::from_utf8_lossy(region.name().as_ref()).into_owned();

        if let Some(ref current_chrom) = chrom {
            if current_chrom != &c {
                return Err(format!("genes are on different chromosomes: {} and {}", current_chrom, c));
            }
        } else {
            chrom = Some(c);
        }

        if let Some(s) = region.interval().start() {
            min_start = min_start.min(usize::from(s));
        }
        if let Some(e) = region.interval().end() {
            max_end = max_end.max(usize::from(e));
        }
    }

    let chrom = chrom.ok_or_else(|| "no regions found for the given genes".to_string())?;

    // Apply offset
    let start = min_start.saturating_sub(offset).max(1);
    let end = max_end.saturating_add(offset);

    let start_pos = Position::try_from(start).map_err(|e| e.to_string())?;
    let end_pos = Position::try_from(end).map_err(|e| e.to_string())?;

    Ok(Region::new(chrom, start_pos..=end_pos))
}
