// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::time::{Instant, Duration};
use std::sync::RwLock;
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};

/// A point-in-time snapshot of memory consumption.
#[derive(Debug, Clone, Copy, Default)]
pub struct MemorySnapshot {
    pub data_bytes: u64,
    pub mapping_bytes: u64,
    pub os_observed_bytes: u64,
    pub total_peak_bytes: u64,
    pub read_count: u64,
    pub header_bytes_uncompressed: u64,
    pub header_bytes_compressed: u64,
    pub payload_bytes_uncompressed: u64,
    pub payload_bytes_compressed: u64,
}

#[derive(Debug, Clone)]
pub struct BamStats {
    pub context: String,
    pub read_count: usize,
    pub segment_count: usize,
    pub data_mb: f64,
    pub mapping_mb: f64,
    pub os_mb: f64,
    pub header_uncompressed_bytes: u64,
    pub header_compressed_bytes: u64,
    pub payload_uncompressed_bytes: u64,
    pub payload_compressed_bytes: u64,
}

pub struct BenchmarkMonitor {
    pub start_time: Instant,
    pub stage_durations: RwLock<HashMap<String, Duration>>,
    pub memory_snapshots: RwLock<HashMap<String, MemorySnapshot>>,
    pub cohort_peaks: RwLock<HashMap<String, u64>>,
    pub current_cohort: RwLock<Option<String>>,
    pub bam_stats: RwLock<Vec<BamStats>>,
    pub total_peak_memory: AtomicU64,
}

impl BenchmarkMonitor {
    pub fn new() -> Self {
        Self {
            start_time: Instant::now(),
            stage_durations: RwLock::new(HashMap::new()),
            memory_snapshots: RwLock::new(HashMap::new()),
            cohort_peaks: RwLock::new(HashMap::new()),
            current_cohort: RwLock::new(None),
            bam_stats: RwLock::new(Vec::new()),
            total_peak_memory: AtomicU64::new(get_os_memory()),
        }
    }

    /// Resets the peak memory tracker to the current baseline usage.
    /// This is useful when switching genes to detect the processing overhead of the new gene.
    pub fn reset_peak_memory(&self) {
        let current_os = get_os_memory();
        self.total_peak_memory.store(current_os, Ordering::SeqCst);
        if let Ok(mut map) = self.cohort_peaks.write() {
            map.clear();
        }
    }

    pub fn set_current_cohort(&self, cohort: Option<String>) {
        if let Ok(mut c) = self.current_cohort.write() {
            *c = cohort;
        }
    }

    pub fn observe_peak(&self) {
        let current_os = get_os_memory();
        self.total_peak_memory.fetch_max(current_os, Ordering::SeqCst);
        
        if let Some(cohort) = self.current_cohort.read().ok().and_then(|c| c.clone()) {
            let key = format!("Cohort: {}", cohort);
            
            // Update dedicated cohort peaks map
            if let Ok(mut map) = self.cohort_peaks.write() {
                let entry = map.entry(cohort.clone()).or_insert(0);
                if current_os > *entry {
                    *entry = current_os;
                }
            }
            
            // Also update the OS component of the active cohort's memory snapshot if it exists
            if let Ok(mut snapshots) = self.memory_snapshots.write() {
                if let Some(s) = snapshots.get_mut(&key) {
                    if current_os > s.os_observed_bytes {
                        s.os_observed_bytes = current_os;
                    }
                    if current_os > s.total_peak_bytes {
                        s.total_peak_bytes = current_os;
                    }
                }
            }
        }
    }

    pub fn record_stage(&self, stage: &str, duration: Duration) {
        if let Ok(mut map) = self.stage_durations.write() {
            map.insert(stage.to_string(), duration);
        }
    }

    pub fn record_memory(&self, context: &str, mut snapshot: MemorySnapshot) {
        if snapshot.os_observed_bytes == 0 {
            snapshot.os_observed_bytes = get_os_memory();
        }
        let peak = snapshot.total_peak_bytes.max(snapshot.os_observed_bytes);
        self.total_peak_memory.fetch_max(peak, Ordering::SeqCst);
        
        let mut key = context.to_string();
        
        // Transform "Building: ..." into "Cohort: ..." if we have an active cohort
        if context.starts_with("Building: ") {
            if let Ok(c_opt) = self.current_cohort.read() {
                if let Some(ref cohort) = *c_opt {
                    key = format!("Cohort: {}", cohort);
                }
            }
        }

        if let Ok(mut map) = self.memory_snapshots.write() {
            // Use entry API to aggregate peaks seen across all snapshots for this key
            let entry = map.entry(key).or_insert(snapshot);
            entry.data_bytes = entry.data_bytes.max(snapshot.data_bytes);
            entry.mapping_bytes = entry.mapping_bytes.max(snapshot.mapping_bytes);
            entry.os_observed_bytes = entry.os_observed_bytes.max(snapshot.os_observed_bytes);
            entry.total_peak_bytes = entry.total_peak_bytes.max(snapshot.total_peak_bytes);
        }

        // Also update dedicated cohort peaks map
        if let Some(cohort) = self.current_cohort.read().ok().and_then(|c| c.clone()) {
            if let Ok(mut map) = self.cohort_peaks.write() {
                let entry = map.entry(cohort).or_insert(0);
                if snapshot.os_observed_bytes > *entry {
                    *entry = snapshot.os_observed_bytes;
                }
            }
        }
    }

    pub fn observe_cohort_peak(&self, cohort: &str) {
        let current_os = get_os_memory();
        self.total_peak_memory.fetch_max(current_os, Ordering::SeqCst);
        if let Ok(mut map) = self.cohort_peaks.write() {
            let entry = map.entry(cohort.to_string()).or_insert(0);
            if current_os > *entry {
                *entry = current_os;
            }
        }
    }

    pub fn record_bam_stats(&self, stats: BamStats) {
        if let Ok(mut list) = self.bam_stats.write() {
            list.push(stats.clone());
        }
        
        // Update cohort-level stats as well
        let key = if let Some(cohort) = self.current_cohort.read().ok().and_then(|c| c.clone()) {
            format!("Cohort: {}", cohort)
        } else {
            "Overall".to_string()
        };
        
        if let Ok(mut snapshots) = self.memory_snapshots.write() {
            let s = snapshots.entry(key).or_default();
            s.read_count += stats.read_count as u64;
            s.header_bytes_uncompressed += stats.header_uncompressed_bytes;
            s.header_bytes_compressed += stats.header_compressed_bytes;
            s.payload_bytes_uncompressed += stats.payload_uncompressed_bytes;
            s.payload_bytes_compressed += stats.payload_compressed_bytes;
        }
    }

    pub fn record_payload_compression(&self, _context_or_path: String, uncompressed: u64, compressed: u64) {
        let key = if let Some(cohort) = self.current_cohort.read().ok().and_then(|c| c.clone()) {
            format!("Cohort: {}", cohort)
        } else {
            "Overall".to_string()
        };

        if let Ok(mut snapshots) = self.memory_snapshots.write() {
            let s = snapshots.entry(key).or_default();
            s.payload_bytes_uncompressed += uncompressed;
            s.payload_bytes_compressed += compressed;
        }
    }

    pub fn write_csv(&self, path: &str) -> std::io::Result<()> {
        use std::fs::File;
        use std::io::Write;
        
        let mut file = File::create(path)?;
        writeln!(file, "Type,Name,Metric,Value,Unit")?;
        
        // Stages
        if let Ok(durations) = self.stage_durations.read() {
            let mut keys: Vec<_> = durations.keys().collect();
            keys.sort();
            for k in keys {
                writeln!(file, "Stage,{},Duration,{:.4},s", k, durations[k].as_secs_f64())?;
            }
        }
        
        // Memory Snapshots
        if let Ok(snapshots) = self.memory_snapshots.read() {
            let mut keys: Vec<_> = snapshots.keys().collect();
            keys.sort();
            for k in keys {
                let s = &snapshots[k];
                writeln!(file, "MemorySnapshot,{},DataMemory,{:.2},MB", k, s.data_bytes as f64 / 1024.0 / 1024.0)?;
                writeln!(file, "MemorySnapshot,{},MappingMemory,{:.2},MB", k, s.mapping_bytes as f64 / 1024.0 / 1024.0)?;
                writeln!(file, "MemorySnapshot,{},OSObserved,{:.2},MB", k, s.os_observed_bytes as f64 / 1024.0 / 1024.0)?;
                writeln!(file, "MemorySnapshot,{},TotalPeak,{:.2},MB", k, s.total_peak_bytes as f64 / 1024.0 / 1024.0)?;
                
                if s.read_count > 0 {
                    writeln!(file, "MemorySnapshot,{},ReadCount,{},count", k, s.read_count)?;
                    let total_mem_bytes = s.data_bytes + s.mapping_bytes;
                    let bytes_per_read = total_mem_bytes as f64 / s.read_count as f64;
                    writeln!(file, "MemorySnapshot,{},MemPerRead,{:.2},bytes", k, bytes_per_read)?;
                    
                    if s.header_bytes_uncompressed > 0 {
                        let rate = 100.0 * (1.0 - (s.header_bytes_compressed as f64 / s.header_bytes_uncompressed as f64));
                        writeln!(file, "MemorySnapshot,{},HeaderCompressionRate,{:.1},%", k, rate)?;
                    }
                    if s.payload_bytes_uncompressed > 0 {
                        let rate = 100.0 * (1.0 - (s.payload_bytes_compressed as f64 / s.payload_bytes_uncompressed as f64));
                        writeln!(file, "MemorySnapshot,{},PayloadCompressionRate,{:.1},%", k, rate)?;
                    }
                }
            }
        }
        
        // BAM Stats
        if let Ok(stats) = self.bam_stats.read() {
            for s in stats.iter() {
                writeln!(file, "RegionStat,{},ReadCount,{},count", s.context, s.read_count)?;
                writeln!(file, "RegionStat,{},SegmentCount,{},count", s.context, s.segment_count)?;
                writeln!(file, "RegionStat,{},DataMem,{:.2},MB", s.context, s.data_mb)?;
                writeln!(file, "RegionStat,{},MappingMem,{:.2},MB", s.context, s.mapping_mb)?;
                writeln!(file, "RegionStat,{},OSMem,{:.2},MB", s.context, s.os_mb)?;
            }
        }

        let total_time = self.start_time.elapsed().as_secs_f64();
        let peak_mem = self.total_peak_memory.load(Ordering::SeqCst) as f64 / 1024.0 / 1024.0 / 1024.0;
        writeln!(file, "Summary,Overall,TotalTime,{:.4},s", total_time)?;
        writeln!(file, "Summary,Overall,PeakMemory,{:.4},GB", peak_mem)?;
        
        Ok(())
    }

    pub fn print_report(&self) {
        println!("\n========================================================");
        println!("                LNC-SEEKER BENCHMARK REPORT             ");
        println!("========================================================");
        
        let durations = self.stage_durations.read().ok();

        println!("\n--- Cohort Performance & Memory Profiles ---");
        if let Ok(snapshots) = self.memory_snapshots.read() {
            let mut keys: Vec<_> = snapshots.keys().collect();
            keys.sort();
            for &k in &keys {
                let s = &snapshots[k];
                let total_mb = s.total_peak_bytes as f64 / 1024.0 / 1024.0;
                let data_mb = s.data_bytes as f64 / 1024.0 / 1024.0;
                let map_mb = s.mapping_bytes as f64 / 1024.0 / 1024.0;
                let os_mb = s.os_observed_bytes as f64 / 1024.0 / 1024.0;
                
                println!("  |-- Context: {}", k);
                
                // Show duration here if it exists for this key
                if let Some(ref d_map) = durations {
                    if let Some(d) = d_map.get(k) {
                        println!("  |   |-- Execution Time:  {:?}", d);
                    } else if k.starts_with("Cohort: ") {
                        if let Some(d) = d_map.get(k.strip_prefix("Cohort: ").unwrap()) {
                            println!("  |   |-- Execution Time:  {:?}", d);
                        }
                    }
                }

                println!("  |   |-- Data Memory:    {:>8.2} MB", data_mb);
                println!("  |   |-- Mapping (Peak): {:>8.2} MB", map_mb);
                println!("  |   |-- OS Observed:    {:>8.2} MB", os_mb);
                println!("  |   |-- Total Peak:     {:>8.2} MB", total_mb);
                
                if s.read_count > 0 {
                    let total_mem_bytes = s.data_bytes + s.mapping_bytes;
                    let bytes_per_read = total_mem_bytes as f64 / s.read_count as f64;
                    println!("  |   |-- Collected Reads: {:>8}", s.read_count);
                    println!("  |   |-- Memory/Read:    {:>8.2} bytes", bytes_per_read);
                    
                    if s.header_bytes_uncompressed > 0 {
                        let h_rate = 100.0 * (1.0 - (s.header_bytes_compressed as f64 / s.header_bytes_uncompressed as f64));
                        println!("  |   |-- Header Comp.:   {:>8.1}%", h_rate);
                    }
                    if s.payload_bytes_uncompressed > 0 {
                        let p_rate = 100.0 * (1.0 - (s.payload_bytes_compressed as f64 / s.payload_bytes_uncompressed as f64));
                        println!("  |   |-- Payload Comp.:  {:>8.1}%", p_rate);
                    }
                }
                println!("  |   L--");
            }
        }

        // Show remaining stages that weren't cohorts
        if let Some(ref d_map) = durations {
            let mut other_stages: Vec<_> = d_map.keys()
                .filter(|k| !k.starts_with("Cohort: ") && !self.memory_snapshots.read().map(|m| m.contains_key(*k)).unwrap_or(false))
                .collect();
            if !other_stages.is_empty() {
                println!("\n--- Other Execution Stages ---");
                other_stages.sort();
                for k in other_stages {
                    println!("  |-- {:<30} : {:?}", k, d_map[k]);
                }
            }
        }

        let peak_gb = self.total_peak_memory.load(Ordering::SeqCst) as f64 / 1024.0 / 1024.0 / 1024.0;
        println!("\nOverall Peak Memory Observed: {:.3} GB", peak_gb);
        println!("Total Elapsed Time: {:?}", self.start_time.elapsed());
        println!("========================================================\n");
    }
}

pub fn get_os_memory() -> u64 {
    #[cfg(target_os = "linux")]
    {
        use std::fs;
        if let Ok(statm) = fs::read_to_string("/proc/self/statm") {
            if let Some(rss_pages_str) = statm.split_whitespace().nth(1) {
                if let Ok(rss_pages) = rss_pages_str.parse::<u64>() {
                    return rss_pages * 4096; // Standard page size
                }
            }
        }
    }

    #[cfg(target_os = "windows")]
    {
        use std::mem;
        use std::os::raw::c_void;

        #[repr(C)]
        #[allow(non_snake_case)]
        struct PROCESS_MEMORY_COUNTERS {
            cb: u32,
            PageFaultCount: u32,
            PeakWorkingSetSize: usize,
            WorkingSetSize: usize,
            QuotaPeakPagedPoolUsage: usize,
            QuotaPagedPoolUsage: usize,
            QuotaPeakNonPagedPoolUsage: usize,
            QuotaNonPagedPoolUsage: usize,
            PagefileUsage: usize,
            PeakPagefileUsage: usize,
        }

        #[link(name = "psapi")]
        unsafe extern "system" {
            fn GetCurrentProcess() -> *mut c_void;
            fn GetProcessMemoryInfo(
                process: *mut c_void,
                counters: *mut PROCESS_MEMORY_COUNTERS,
                size: u32,
            ) -> i32;
        }

        let mut counters: PROCESS_MEMORY_COUNTERS = unsafe { mem::zeroed() };
        counters.cb = mem::size_of::<PROCESS_MEMORY_COUNTERS>() as u32;
        let res = unsafe { GetProcessMemoryInfo(GetCurrentProcess(), &mut counters, counters.cb) };
        if res != 0 {
            return counters.WorkingSetSize as u64;
        }
    }

    0
}
