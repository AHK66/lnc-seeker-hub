// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::path::{Path, PathBuf};
use std::collections::HashSet;
pub use lnc_seeker_lib::pipeline_config::{AnyResult, Cohort};
use lnc_seeker_lib::pipeline_config::parse_recursive;
use std::cell::RefCell;

#[derive(Debug, Clone)]
pub struct Config {
    pub gtf_path: String,
    pub assembly_report: String,
    pub genome_path: Option<String>,
    pub cohorts: Vec<Cohort>,
    pub output_dir: Option<String>,
    pub path_prefix: Option<String>,
    pub annotate_endpoint_tag: bool,
    pub skip_secondary: bool,
    pub write_bam: bool,
    pub extra_regions: Vec<String>,
    pub gene_regions: Vec<(Vec<String>, usize)>,
    pub compress_headers: bool,
    pub compression_algorithm: String,
    pub compression_use_substitutes: bool,
    pub use_delta_incremental_collection: bool,
    pub stream: bool,
}

pub fn parse_config(config_path: &str, extra_includes: &[String]) -> AnyResult<Config> {
    let context = RefCell::new(ParserContext {
        gtf_path: String::new(),
        assembly_report: String::new(),
        genome_path: None,
        annotate_endpoint_tag: false,
        skip_secondary: true,
        write_bam: false,
        extra_regions: Vec::new(),
        gene_regions: Vec::new(),
        compress_headers: false,
        compression_algorithm: "huffman".to_string(),
        compression_use_substitutes: true,
        use_delta_incremental_collection: false,
        stream: true,
        cohorts: Vec::new(),
        current_cohort_name: "default".to_string(),
        current_bam_paths: Vec::new(),
        current_output_bam: None,
        current_tissue: None,
        current_status: None,
        output_dir: None,
        path_prefix: None,
        visited_files: HashSet::new(),
    });

    let abs_path = std::fs::canonicalize(config_path)?;
    
    let mut visited_files = HashSet::new();

    {
        let mut wrap_handler = |key: &str, value: &str, path: &Path| -> AnyResult<()> {
            context.borrow_mut().handle_line(key, value, path)
        };
        parse_recursive(&abs_path, &mut visited_files, 0, &mut wrap_handler)?;

        let base_abs_path = abs_path.parent().unwrap_or(Path::new("."));
        for include_file in extra_includes {
            let abs_include_path = base_abs_path.join(include_file);
            parse_recursive(&abs_include_path, &mut visited_files, 0, &mut wrap_handler)?;
        }
    }

    let mut context = context.into_inner();
    context.visited_files = visited_files;
    context.push_current_cohort()?;

    if context.gtf_path.is_empty() { return Err("config missing `gtf=` entry".into()); }
    if context.assembly_report.is_empty() { return Err("config missing `assembly_report=` entry".into()); }
    if context.cohorts.is_empty() { return Err("config must contain at least one cohort".into()); }

    Ok(Config {
        gtf_path: context.gtf_path,
        assembly_report: context.assembly_report,
        genome_path: context.genome_path,
        cohorts: context.cohorts,
        output_dir: context.output_dir,
        path_prefix: context.path_prefix,
        annotate_endpoint_tag: context.annotate_endpoint_tag,
        skip_secondary: context.skip_secondary,
        write_bam: context.write_bam,
        extra_regions: context.extra_regions,
        gene_regions: context.gene_regions,
        compress_headers: context.compress_headers,
        compression_algorithm: context.compression_algorithm,
        compression_use_substitutes: context.compression_use_substitutes,
        use_delta_incremental_collection: context.use_delta_incremental_collection,
        stream: context.stream,
    })
}


struct ParserContext {
    gtf_path: String,
    assembly_report: String,
    genome_path: Option<String>,
    annotate_endpoint_tag: bool,
    skip_secondary: bool,
    write_bam: bool,
    extra_regions: Vec<String>,
    gene_regions: Vec<(Vec<String>, usize)>,
    compress_headers: bool,
    compression_algorithm: String,
    compression_use_substitutes: bool,
    use_delta_incremental_collection: bool,
    stream: bool,
    cohorts: Vec<Cohort>,
    current_cohort_name: String,
    current_bam_paths: Vec<String>,
    current_output_bam: Option<String>,
    current_tissue: Option<String>,
    current_status: Option<String>,
    output_dir: Option<String>,
    path_prefix: Option<String>,
    visited_files: HashSet<PathBuf>,
}

impl ParserContext {
    fn handle_line(&mut self, key: &str, value: &str, _path: &Path) -> AnyResult<()> {
        match key {
            "gtf" => self.gtf_path = value.to_string(),
            "assembly_report" => self.assembly_report = value.to_string(),
            "genome" => self.genome_path = Some(value.to_string()),
            "annotate_endpoint_tag" => self.annotate_endpoint_tag = value.parse().unwrap_or(false),
            "skip_secondary" => self.skip_secondary = value.parse().unwrap_or(true),
            "write_bam" => self.write_bam = value.parse().unwrap_or(false),
            "output_dir" => self.output_dir = Some(value.to_string()),
            "path_prefix" => self.path_prefix = Some(value.to_string()),
            "bam" => self.current_bam_paths.push(value.to_string()),
            "output_bam" => self.current_output_bam = Some(value.to_string()),
            "tissue" => self.current_tissue = Some(value.to_string()),
            "status" => self.current_status = Some(value.to_string()),
            "region" => self.extra_regions.push(value.to_string()),
            "gene_region" => {
                let (names_str, offset_str) = value.split_once(':').unwrap_or((value, "0"));
                let names: Vec<String> = names_str.split(',').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect();
                let offset = offset_str.parse::<usize>().unwrap_or(0);
                if !names.is_empty() { self.gene_regions.push((names, offset)); }
            },
            "compress_headers" => self.compress_headers = value.parse().unwrap_or(false),
            "compression_algorithm" => self.compression_algorithm = value.to_string(),
            "compression_use_substitutes" => self.compression_use_substitutes = value.parse().unwrap_or(true),
            "use_delta_incremental_collection" => self.use_delta_incremental_collection = value.parse().unwrap_or(false),
            "stream" => self.stream = value.parse().unwrap_or(true),
            "section" => {
                if value.starts_with("cohort:") {
                    self.push_current_cohort()?;
                    self.current_cohort_name = value.strip_prefix("cohort:").unwrap().to_string();
                }
            },
            _ => { /* Ignore unknown keys */ }
        }
        Ok(())
    }

    fn push_current_cohort(&mut self) -> AnyResult<()> {
        if !self.current_bam_paths.is_empty() {
            let bam_paths: Vec<String> = self.current_bam_paths.drain(..).collect();
            let num_samples = bam_paths.len();
            self.cohorts.push(Cohort {
                name: self.current_cohort_name.clone(),
                bam_paths,
                output_bam: self.current_output_bam.take(),
                tissue: self.current_tissue.take(),
                status: self.current_status.take(),
                num_samples,
            });
        }
        Ok(())
    }
}
