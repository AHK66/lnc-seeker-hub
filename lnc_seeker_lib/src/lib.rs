// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
// PyO3 0.21 macros generate code that triggers unsafe_op_in_unsafe_fn warnings in Rust 2024.
// This is allowed until PyO3 adds explicit unsafe blocks to its generated expansion.
#![allow(unsafe_op_in_unsafe_fn)]

use std::sync::Arc;
use std::time::Instant;
use pyo3::prelude::*;

pub mod assembly;
pub mod utils;
pub mod progress;
pub mod config;
pub mod models;
pub mod genome;
pub mod coverage;
pub mod annotations;
pub mod analysis;
pub mod processing;
pub mod reads_manager;
pub mod compression;
pub mod regions;
pub mod gtf;
pub mod pipeline_config;
#[path = "../../benchmarking/rust/benchmark_core.rs"]
pub mod benchmarking;

pub use benchmarking::{BenchmarkMonitor, MemorySnapshot};

#[cfg(test)]
pub mod tests;

pub use assembly::AssemblyReport;
pub use utils::{FeatureAttributes, parse_attributes, intersects};
pub use progress::{ProgressData, SessionProgress, get_session_progress};
pub use config::Config;
pub use models::{SampleResult, AnalysisResult, JunctionSpan, ReadInfo, Annotation};
pub use analysis::{run_analysis, cleanup_csv_files};
pub use annotations::optimize_gtf_to_lba;

#[pyfunction]
fn load_config_py(path: String) -> PyResult<String> {
    std::fs::read_to_string(path).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
}

#[pyfunction]
fn run_analysis_py(py: Python<'_>, config_json: String, progress: PyRef<SessionProgress>) -> PyResult<String> {
    let t_begin = Instant::now();
    let config: Config = serde_json::from_str(&config_json).map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
    println!("[BENCHMARK] Rust Binding: Deserialization took {:?}", t_begin.elapsed());

    let progress_data = Arc::clone(&progress.data);
    let t_analysis_start = Instant::now();
    let result = py.allow_threads(move || {
        run_analysis(&config, progress_data)
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    println!("[BENCHMARK] Rust Binding: Core analysis took {:?}", t_analysis_start.elapsed());
    
    let t_ser_start = Instant::now();
    let res = serde_json::to_string(&result).map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
    println!("[BENCHMARK] Rust Binding: Serialization took {:?}", t_ser_start.elapsed());
    println!("[BENCHMARK] Rust Binding: Total wrapper time {:?}", t_begin.elapsed());
    Ok(res)
}

#[pyfunction]
fn filter_outliers_py(sample_json: String) -> PyResult<String> {
    let mut sample: SampleResult = serde_json::from_str(&sample_json).map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
    sample.filter_outliers();
    serde_json::to_string(&sample).map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}

#[pyfunction]
fn get_junction_reads_py(py: Python<'_>, bam_path: String, reference: String, start_target: usize, end_target: usize, min_mq: u8, max_reads: usize, genome_path: Option<String>, filter_clean: Option<bool>) -> PyResult<String> {
    let filter_clean_val = filter_clean.unwrap_or(false);
    let result = py.allow_threads(move || {
        processing::get_junction_reads(&bam_path, &reference, start_target, end_target, min_mq, max_reads, &genome_path, filter_clean_val)
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    serde_json::to_string(&result).map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}

#[pyfunction]
fn get_junction_reads_batch_py(
    py: Python<'_>,
    bam_path: String,
    reference: String,
    junctions: Vec<(usize, usize)>,
    min_mq: u8,
    max_reads: usize,
    genome_path: Option<String>,
    filter_clean: Option<bool>,
) -> PyResult<String> {
    let filter_clean_val = filter_clean.unwrap_or(false);
    let t_start = Instant::now();
    let result = py
        .allow_threads(move || {
            processing::get_junction_reads_batch(
                &bam_path,
                &reference,
                &junctions,
                min_mq,
                max_reads,
                &genome_path,
                filter_clean_val,
            )
        })
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    let read_provider = crate::reads_manager::get_read_provider();
    read_provider.benchmark.record_stage("Batch Read Fetch", t_start.elapsed());
    read_provider.benchmark.observe_peak();

    serde_json::to_string(&result).map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}


#[pyfunction]
fn downsample_coverage_py(
    positions: Vec<usize>,
    depths_bg: Vec<f64>,
    depths_fg: Vec<f64>,
    target_points: usize,
    forced_positions: Vec<usize>,
) -> (Vec<usize>, Vec<f64>, Vec<f64>) {
    let t_start = Instant::now();
    let res = processing::downsample_coverage(positions, depths_bg, depths_fg, target_points, forced_positions);
    
    let read_provider = crate::reads_manager::get_read_provider();
    read_provider.benchmark.record_stage("Downsample Coverage", t_start.elapsed());
    
    res
}

#[pyfunction]
fn clear_all_caches_py() {
    crate::reads_manager::get_read_provider().clear_all_caches();
}

#[pyfunction]
fn get_cache_core_size_py() -> f64 {
    crate::reads_manager::get_read_provider().get_core_cache_usage_mb()
}

#[pyfunction]
fn get_cache_annotation_size_py() -> f64 {
    crate::reads_manager::get_read_provider().get_annotation_cache_usage_mb()
}

#[pymodule]
fn lnc_seeker(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<SessionProgress>()?;
    m.add_function(wrap_pyfunction!(load_config_py, m)?)?;
    m.add_function(wrap_pyfunction!(run_analysis_py, m)?)?;
    m.add_function(wrap_pyfunction!(filter_outliers_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_junction_reads_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_junction_reads_batch_py, m)?)?;
    m.add_function(wrap_pyfunction!(downsample_coverage_py, m)?)?;
    m.add_function(wrap_pyfunction!(clear_all_caches_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_cache_core_size_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_cache_annotation_size_py, m)?)?;
    Ok(())
}
