// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::path::{Path, PathBuf};
use std::collections::HashSet;
use std::fs;
use std::io::{BufRead, BufReader};

pub type AnyResult<T> = Result<T, Box<dyn std::error::Error>>;

/// A collection of BAM files sharing metadata.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Cohort {
    pub name: String,
    pub bam_paths: Vec<String>,
    pub output_bam: Option<String>,
    pub tissue: Option<String>,
    pub status: Option<String>,
    pub num_samples: usize,
}

/// Generic pipeline configuration parser infrastructure.
pub struct BaseParserContext {
    pub visited_files: HashSet<PathBuf>,
}

impl BaseParserContext {
    pub fn new() -> Self {
        Self {
            visited_files: HashSet::new(),
        }
    }
}

/// Recursively parses a key-value style configuration file.
pub fn parse_recursive<F>(
    config_path: &Path,
    visited_files: &mut HashSet<PathBuf>,
    depth: usize,
    mut line_handler: F,
) -> AnyResult<()>
where
    F: FnMut(&str, &str, &Path) -> AnyResult<()>,
{
    parse_recursive_impl(config_path, visited_files, depth, &mut line_handler)
}

fn parse_recursive_impl(
    config_path: &Path,
    visited_files: &mut HashSet<PathBuf>,
    depth: usize,
    line_handler: &mut dyn FnMut(&str, &str, &Path) -> AnyResult<()>,
) -> AnyResult<()> {
    if depth > 10 {
        return Err("Too many recursive includes".into());
    }

    let config_path = fs::canonicalize(config_path)?;
    if !visited_files.insert(config_path.clone()) {
        return Ok(());
    }

    let file = fs::File::open(&config_path)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if let Some(include_path_str) = trimmed.strip_prefix("include=") {
            let base_dir = config_path.parent().unwrap_or(Path::new("."));
            let include_path = base_dir.join(include_path_str);
            parse_recursive_impl(&include_path, visited_files, depth + 1, line_handler)?;
            continue;
        }

        if let Some((key, value)) = trimmed.split_once('=') {
            line_handler(key.trim(), value.trim(), &config_path)?;
        } else if trimmed.starts_with('[') && trimmed.ends_with(']') {
            line_handler("section", &trimmed[1..trimmed.len() - 1], &config_path)?;
        }
    }

    Ok(())
}
