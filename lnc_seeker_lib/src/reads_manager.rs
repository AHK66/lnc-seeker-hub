// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::collections::HashMap;
use std::sync::{RwLock, OnceLock, Arc};
use std::time::Instant;
use std::fs::File;
use std::io::{Read, Write, BufWriter};
use std::path::PathBuf;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use serde::{Serialize, Deserialize};
use crate::models::{ReadInfo, ReadSegment, LbaFile};
use crate::utils::normalize_path;
use crate::compression::{CompressionMode, EditOp, compute_name_diff};
use crate::{BenchmarkMonitor, MemorySnapshot};

const CACHE_MAGIC_V1: &[u8; 4] = b"LNC1";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CacheStatus {
    Empty,
    Caching,
    Complete,
    ExceededLimit,
}

/// A highly compact representation of a BAM alignment record.
/// Takes only ~24 bytes plus amortized segment/name storage.
#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct CompactRead {
    pub name_id: u32,
    pub ref_id: u16,
    pub mq: u8,
    pub flags: u8, // bit 0: is_paired, bit 1: is_second_mate, bits 2-3: strand (0 unknown, 1 plus, 2 minus)
    pub start: u32,
    pub end: u32,
    pub seg_offset: u32,
    pub seg_count: u8,
}

#[inline]
fn encode_strand_bits(strand: i8) -> u8 {
    match strand {
        1 => 0x04,
        -1 => 0x08,
        _ => 0x00,
    }
}

#[inline]
fn decode_strand_bits(flags: u8) -> i8 {
    match (flags >> 2) & 0x03 {
        1 => 1,
        2 => -1,
        _ => 0,
    }
}

#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct CompactSegment {
    pub start: u32,
    pub end: u32,
}

/// Incremental delta-encoded name storage to avoid memory peaks during scanning.
#[derive(Clone, Serialize, Deserialize)]
pub struct DeltaNameStore {
    pub symbols: Vec<u16>,
    pub num_names: usize,
    pub use_substitutes: bool,
    #[serde(skip)]
    pub last_name: Vec<u8>,
}

impl DeltaNameStore {
    pub fn new(use_substitutes: bool) -> Self {
        Self {
            symbols: Vec::with_capacity(1_000_000),
            num_names: 0,
            use_substitutes,
            last_name: Vec::new(),
        }
    }

    pub fn add_name(&mut self, name: &str) {
        let current = name.as_bytes();
        let ops = compute_name_diff(&self.last_name, current, self.use_substitutes);
        let end_symbol = if self.use_substitutes { 1024 } else { 768 };

        for op in ops {
            match op {
                EditOp::Match(len) => self.symbols.push(256 + len as u16),
                EditOp::Delete(len) => self.symbols.push(512 + len as u16),
                EditOp::Substitute(data) => {
                    self.symbols.push(768 + data.len() as u16);
                    for &b in &data { self.symbols.push(b as u16); }
                }
                EditOp::Insert(data) => {
                    for &b in &data { self.symbols.push(b as u16); }
                }
            }
        }
        self.symbols.push(end_symbol);
        self.last_name = current.to_vec();
        self.num_names += 1;
    }

    pub fn to_names(&self) -> Vec<String> {
        if self.num_names == 0 { return Vec::new(); }
        crate::compression::symbols_to_names(&self.symbols, self.num_names, self.use_substitutes)
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct BamCache {
    pub reads: Vec<CompactRead>,
    pub segments: Vec<CompactSegment>,
    pub segment_tags: Vec<u8>,
    pub names: Vec<String>,
    #[serde(skip)]
    pub delta_names: Option<DeltaNameStore>,
    pub ref_names: Vec<String>,
    pub ref_lengths: Vec<usize>,
    pub status: CacheStatus,
    #[serde(skip)]
    pub compressed_headers: Option<crate::compression::CompressedHeaders>,
    #[serde(skip)]
    pub size_bytes: usize,
    
    // Telemetry for compression
    #[serde(skip)]
    pub zstd_payload_uncompressed: u64,
    #[serde(skip)]
    pub zstd_payload_compressed: u64,
}

impl BamCache {
    pub fn get_uncompressed_header_bytes(&self) -> u64 {
        if let Some(ref d) = self.delta_names {
            // Estimate based on symbols (rough)
            (d.symbols.len() * 2) as u64
        } else {
            self.names.iter().map(|s| s.len() as u64).sum()
        }
    }

    pub fn get_compressed_header_bytes(&self) -> u64 {
        self.compressed_headers.as_ref().map(|c| c.data.len() as u64).unwrap_or(0)
    }

    pub fn get_payload_bytes(&self) -> u64 {
        (self.reads.len() * std::mem::size_of::<CompactRead>() + 
         self.segments.len() * std::mem::size_of::<CompactSegment>() +
         self.segment_tags.len()) as u64
    }

    /// Calculates the sum of all segment lengths in the cache.
    pub fn total_coverage_bases(&self) -> u64 {
        self.segments.iter().map(|s| (s.end - s.start) as u64).sum()
    }

    pub fn recalculate_size_bytes(&mut self) {
        let mem_reads = self.reads.capacity() * std::mem::size_of::<CompactRead>();
        let mem_segs = self.segments.capacity() * std::mem::size_of::<CompactSegment>();
        let mem_tags = self.segment_tags.capacity();
        let mem_names: usize = self.names.capacity() * std::mem::size_of::<String>() + 
                               self.names.iter().map(|s| s.capacity()).sum::<usize>();
        let mem_delta = self.delta_names.as_ref().map(|d| {
            d.symbols.capacity() * std::mem::size_of::<u16>() + d.last_name.capacity()
        }).unwrap_or(0);
        let mem_ref_names: usize = self.ref_names.capacity() * std::mem::size_of::<String>() + 
                                   self.ref_names.iter().map(|s| s.capacity()).sum::<usize>();
        let mem_ref_lengths = self.ref_lengths.capacity() * std::mem::size_of::<usize>();
        
        let mem_compressed = if let Some(ref c) = self.compressed_headers {
            std::mem::size_of::<crate::compression::CompressedHeaders>() + 
            c.data.capacity() + 
            c.symbol_lengths.as_ref().map(|s| s.capacity()).unwrap_or(0)
        } else {
            0
        };

        self.size_bytes = mem_reads + mem_segs + mem_tags + mem_names + mem_delta + mem_ref_names + mem_ref_lengths + mem_compressed;
    }

    fn materialized_names(&self) -> Vec<String> {
        if !self.names.is_empty() {
            return self.names.clone();
        }
        if let Some(ref delta) = self.delta_names {
            return delta.to_names();
        }
        if let Some(ref compressed) = self.compressed_headers {
            return crate::compression::decompress_header_block(compressed);
        }
        Vec::new()
    }

    pub fn compress_headers(&mut self) {
        if self.names.is_empty() && self.delta_names.is_none() {
            return;
        }
        if self.compressed_headers.is_some() {
            return;
        }

        let before_mb = self.size_bytes as f64 / 1024.0 / 1024.0;
        let provider = get_read_provider();
        let mode = provider.compression_mode.read().unwrap_or_else(|e| e.into_inner())
            .unwrap_or(CompressionMode::Huffman);
        let use_substitutes = *provider.use_substitutes.read().unwrap_or_else(|e| e.into_inner());

        let compressed = if let Some(delta) = self.delta_names.take() {
            // Re-encode symbols if mode is Huffman or Zstd to ensure optimal distribution
            // and include literal bytes properly. Note: delta symbols are already close to None mode.
            let names = delta.to_names();
            let c = crate::compression::compress_header_block(&names, mode, use_substitutes);
            c
        } else {
            crate::compression::compress_header_block(&self.names, mode, use_substitutes)
        };

        self.compressed_headers = Some(compressed);
        self.names.clear();
        self.names.shrink_to_fit();
        self.recalculate_size_bytes();

        let after_mb = self.size_bytes as f64 / 1024.0 / 1024.0;
        println!("Rust: Post-processing compression for names complete. Memory: {:.2} MB -> {:.2} MB ({:.1}% reduction)", 
                 before_mb, after_mb, (1.0 - after_mb / before_mb.max(0.001)) * 100.0);
    }
}

pub struct IncrementalCacheBuilder {
    name_to_id: HashMap<String, u32>,
    name_hash_to_id: HashMap<u128, u32>,
    ref_to_id: HashMap<String, u16>,
    pub cache: BamCache,
    pub bam_path: String,
}

fn hash_name_128(name: &str) -> u128 {
    let mut h1 = DefaultHasher::new();
    name.hash(&mut h1);
    let v1 = h1.finish();
    let mut h2 = DefaultHasher::new();
    "delta_salt".hash(&mut h2);
    name.hash(&mut h2);
    let v2 = h2.finish();
    ((v1 as u128) << 64) | (v2 as u128)
}

impl IncrementalCacheBuilder {
    pub fn new(bam_path: &str, ref_names: Vec<String>, ref_lengths: Vec<usize>) -> Self {
        let ref_to_id = ref_names.iter().enumerate()
            .map(|(i, n)| (n.clone(), i as u16)).collect();
        
        let provider = get_read_provider();
        let use_delta = *provider.use_delta_incremental.read().unwrap_or_else(|e| e.into_inner());
        let use_substitutes = *provider.use_substitutes.read().unwrap_or_else(|e| e.into_inner());

        let delta_names = if use_delta {
            Some(DeltaNameStore::new(use_substitutes))
        } else {
            None
        };

        Self {
            name_to_id: if use_delta { HashMap::new() } else { HashMap::with_capacity(100_000) },
            name_hash_to_id: if use_delta { HashMap::with_capacity(100_000) } else { HashMap::new() },
            ref_to_id,
            bam_path: bam_path.to_string(),
            cache: BamCache {
                reads: Vec::new(),
                segments: Vec::new(),
                segment_tags: Vec::new(),
                names: Vec::new(),
                delta_names,
                ref_names,
                ref_lengths,
                status: CacheStatus::Caching,
                compressed_headers: None,
                size_bytes: 0,
                zstd_payload_uncompressed: 0,
                zstd_payload_compressed: 0,
            },
        }
    }

    pub fn get_mapping_memory_estimate(&self) -> usize {
        let mem_id = self.name_to_id.capacity() * (std::mem::size_of::<String>() + std::mem::size_of::<u32>()) +
                     self.name_to_id.keys().map(|s| s.capacity()).sum::<usize>();
        let mem_hash = self.name_hash_to_id.capacity() * (std::mem::size_of::<u128>() + std::mem::size_of::<u32>());
        let mem_ref = self.ref_to_id.capacity() * (std::mem::size_of::<String>() + std::mem::size_of::<u16>()) +
                      self.ref_to_id.keys().map(|s| s.capacity()).sum::<usize>();
        mem_id + mem_hash + mem_ref
    }

    pub fn record_snapshot(&mut self) {
        self.cache.recalculate_size_bytes();
        let data_bytes = self.cache.size_bytes;
        let mapping_bytes = self.get_mapping_memory_estimate();
        
        get_read_provider().benchmark.record_memory(
            &format!("Building: {}", self.bam_path),
            MemorySnapshot {
                data_bytes: data_bytes as u64,
                mapping_bytes: mapping_bytes as u64,
                os_observed_bytes: 0,
                total_peak_bytes: (data_bytes + mapping_bytes) as u64,
                ..Default::default()
            }
        );
    }

    pub fn add_read(&mut self, read: ReadInfo) {
        let is_new_name;
        let name_id = if let Some(ref mut delta) = self.cache.delta_names {
            let hash = hash_name_128(&read.name);
            if let Some(&id) = self.name_hash_to_id.get(&hash) {
                is_new_name = false;
                id
            } else {
                let id = delta.num_names as u32;
                delta.add_name(&read.name);
                self.name_hash_to_id.insert(hash, id);
                is_new_name = true;
                id
            }
        } else {
            if let Some(&id) = self.name_to_id.get(&read.name) {
                is_new_name = false;
                id
            } else {
                let id = self.cache.names.len() as u32;
                let name_entry = read.name.clone();
                self.cache.names.push(read.name);
                self.name_to_id.insert(name_entry, id);
                is_new_name = true;
                id
            }
        };

        // Record snapshot every 25,000 new names to track peak accurately
        if is_new_name && (self.name_hash_to_id.len() + self.name_to_id.len()) % 25_000 == 0 {
            self.record_snapshot();
        }

        let ref_id = *self.ref_to_id.get(&read.reference).unwrap_or(&0);
        
        let mut flags = 0u8;
        flags |= 0x01; // paired
        if read.segments.first().map(|s| s.is_mate).unwrap_or(false) {
            flags |= 0x02;
        }
        flags |= encode_strand_bits(read.strand);

        let seg_offset = self.cache.segments.len() as u32;
        let seg_count = read.segments.len() as u8;
        
        for s in read.segments {
            self.cache.segments.push(CompactSegment { start: s.start as u32, end: s.end as u32 });
            
            let mut tag = 0u8;
            if s.is_followed_by_deletion { tag |= 0x01; }
            tag |= (s.mismatches & 0x0F) << 1;
            tag |= (s.insertions & 0x07) << 5;
            self.cache.segment_tags.push(tag);
        }

        self.cache.reads.push(CompactRead {
            name_id,
            ref_id,
            mq: read.mapping_quality,
            flags,
            start: read.start as u32,
            end: read.end as u32,
            seg_offset,
            seg_count,
        });
    }

    pub fn finalize(mut self, status: CacheStatus) -> BamCache {
        self.record_snapshot(); // Final snapshot before dropping builder maps
        self.cache.status = status;
        self.cache.recalculate_size_bytes();
        self.cache
    }
}

pub struct ReadProvider {
    pub caches: RwLock<HashMap<String, BamCache>>,
    pub lba_caches: RwLock<HashMap<String, Arc<LbaFile>>>,
    pub max_memory_mb: RwLock<f64>,
    pub compression_mode: RwLock<Option<CompressionMode>>,
    pub use_substitutes: RwLock<bool>,
    pub use_delta_incremental: RwLock<bool>,
    pub benchmark: Arc<BenchmarkMonitor>,
}

pub static READ_PROVIDER: OnceLock<ReadProvider> = OnceLock::new();

pub fn get_read_provider() -> &'static ReadProvider {
    READ_PROVIDER.get_or_init(|| ReadProvider {
        caches: RwLock::new(HashMap::new()),
        lba_caches: RwLock::new(HashMap::new()),
        max_memory_mb: RwLock::new(2048.0),
        compression_mode: RwLock::new(None),
        use_substitutes: RwLock::new(true),
        use_delta_incremental: RwLock::new(false),
        benchmark: Arc::new(BenchmarkMonitor::new()),
    })
}

impl ReadProvider {
    pub fn set_compress_headers(&self, mode: Option<CompressionMode>, use_substitutes: bool) {
        if let Ok(mut m) = self.compression_mode.write() {
            *m = mode;
        }
        if let Ok(mut s) = self.use_substitutes.write() {
            *s = use_substitutes;
        }
    }

    pub fn set_delta_incremental(&self, enabled: bool) {
        if let Ok(mut d) = self.use_delta_incremental.write() {
            *d = enabled;
        }
    }

    pub fn set_max_memory_mb(&self, mb: f64) {
        if let Ok(mut max) = self.max_memory_mb.write() {
            *max = mb;
        }
    }

    pub fn get_core_cache_usage_mb(&self) -> f64 {
        let bam_bytes: usize = self.caches.read().unwrap_or_else(|e| e.into_inner())
            .values().map(|c| c.size_bytes).sum();
        bam_bytes as f64 / 1024.0 / 1024.0
    }

    pub fn get_annotation_cache_usage_mb(&self) -> f64 {
        let lba_bytes: usize = self.lba_caches.read().unwrap_or_else(|e| e.into_inner())
            .values().map(|l| l.estimate_size_bytes()).sum();
        lba_bytes as f64 / 1024.0 / 1024.0
    }

    pub fn get_total_memory_usage_mb(&self) -> f64 {
        self.get_core_cache_usage_mb() + self.get_annotation_cache_usage_mb()
    }

    /// Clears all LNC1 and LBA caches and resets the peak memory tracker.
    pub fn clear_all_caches(&self) {
        if let Ok(mut caches) = self.caches.write() {
            caches.clear();
        }
        if let Ok(mut lba_caches) = self.lba_caches.write() {
            lba_caches.clear();
        }
        self.benchmark.reset_peak_memory();
    }

    pub fn is_in_ram(&self, bam_path: &str) -> bool {
        let normalized = normalize_path(bam_path);
        if let Ok(caches) = self.caches.read() {
            caches.contains_key(&normalized)
        } else {
            false
        }
    }

    pub fn retain_only_selected(&self, selected_paths: &[String]) {
        let normalized_selected: std::collections::HashSet<String> = selected_paths.iter()
            .map(|p| normalize_path(p))
            .collect();
        
        if let Ok(mut caches) = self.caches.write() {
            let before = caches.len();
            caches.retain(|path, _| normalized_selected.contains(path));
            let after = caches.len();
            if before != after {
                println!("Rust Manager: Purged {} unselected caches from RAM ({} remaining)", before - after, after);
            }
        }
        self.print_cache_layout();
    }

    pub fn print_cache_layout(&self) {
        if let Ok(caches) = self.caches.read() {
            println!("Rust Read Cache Layout:");
            // Collect and sort keys for consistent output
            let mut keys: Vec<_> = caches.keys().collect();
            keys.sort();

            for path in keys {
                let cache = &caches[path];
                let mb = cache.size_bytes as f64 / 1024.0 / 1024.0;
                let status_str = match cache.status {
                    CacheStatus::Complete => "Complete",
                    CacheStatus::Caching => {
                        if cache.delta_names.is_some() { "Caching (Delta)..." } else { "Caching..." }
                    },
                    CacheStatus::ExceededLimit => "ExceededLimit",
                    CacheStatus::Empty => "Empty",
                };
                if cache.size_bytes > 0 || cache.status == CacheStatus::Caching || cache.status == CacheStatus::ExceededLimit {
                   println!("  |-- {}: {:.2} MB ({})", path, mb, status_str);
                }
            }
        }

        if let Ok(lba_caches) = self.lba_caches.read() {
            if !lba_caches.is_empty() {
                println!("Rust LBA Cache Layout:");
                for (path, lba) in lba_caches.iter() {
                    let mb = lba.estimate_size_bytes() as f64 / 1024.0 / 1024.0;
                    println!("  |-- {}: {:.2} MB ({} records)", path, mb, lba.records.len());
                }
            }
        }

        let current_mb = self.get_total_memory_usage_mb();
        let limit = *self.max_memory_mb.read().unwrap_or_else(|e| e.into_inner());
        println!("  L-- Total: {:.2} MB / {:.1} MB", current_mb, limit);
    }

    pub fn get_or_load_lba(&self, path: &str) -> std::io::Result<Arc<LbaFile>> {
        let normalized = normalize_path(path);
        
        // 1. Try RAM
        if let Ok(caches) = self.lba_caches.read() {
            if let Some(lba) = caches.get(&normalized) {
                return Ok(Arc::clone(lba));
            }
        }
        
        // 2. Load Disk
        let t_start = Instant::now();
        let file = File::open(path)?;
        let mut reader = std::io::BufReader::new(file);
        let lba: LbaFile = bincode::deserialize_from(&mut reader)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        
        let arc_lba = Arc::new(lba);
        
        // 3. Store in RAM
        if let Ok(mut caches) = self.lba_caches.write() {
            caches.insert(normalized, Arc::clone(&arc_lba));
        }
        
        println!("[BENCHMARK] Rust Manager: LBA disk-to-RAM load took {:?}", t_start.elapsed());
        self.print_cache_layout();
        
        Ok(arc_lba)
    }

    fn get_cache_path(&self, bam_path: &str) -> PathBuf {
        let normalized = normalize_path(bam_path);
        if normalized.ends_with(".lnc_cache.bin") {
             return PathBuf::from(normalized);
        }
        let mut path = PathBuf::from(normalized);
        path.set_extension("lnc_cache.bin");
        path
    }

    pub fn persist_to_disk(&self, bam_path: &str) -> std::io::Result<()> {
        let normalized = normalize_path(bam_path);
        let caches = self.caches.read().map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        
        if let Some(cache) = caches.get(&normalized) {
            // Only persist complete caches
            if cache.status != CacheStatus::Complete {
                return Ok(());
            }

            let cache_path = self.get_cache_path(bam_path);
            let mut file = File::create(cache_path)?;
            
            // New multi-block format starting with magic V1
            file.write_all(CACHE_MAGIC_V1)?;

            let mut writer = BufWriter::new(file);

            // Block 1: Raw (already compressed) Huffman data
            bincode::serialize_into(&mut writer, &cache.compressed_headers)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

            // Block 2: Zstd compressed remainder (excluding the bitstream)
            // Measure raw size first
            let raw_zstd_data = bincode::serialize(&cache)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            let uncompressed_len = raw_zstd_data.len() as u64;

            let mut inner_writer = Vec::new();
            {
                let mut encoder = zstd::stream::Encoder::new(&mut inner_writer, 3)?;
                encoder.write_all(&raw_zstd_data)?;
                encoder.finish()?;
            }
            let compressed_len = inner_writer.len() as u64;
            
            writer.write_all(&inner_writer)?;
            writer.flush()?;

            // Update cache telemetry (requires mutable access which we usually have during write)
            // Since we're in a read lock for the hashmap, we might need to update a side-channel 
            // but we can at least record it for the benchmark now.
            let provider = get_read_provider();
            provider.benchmark.record_payload_compression(normalized.clone(), uncompressed_len, compressed_len);
        }
        Ok(())
    }

    pub fn load_from_disk(&self, bam_path: &str) -> std::io::Result<bool> {
        let t_start = Instant::now();
        let cache_path = self.get_cache_path(bam_path);
        if !cache_path.exists() {
            return Ok(false);
        }

        let mut file = File::open(&cache_path)?;
        let file_size = file.metadata()?.len();

        // Read entire file into memory for maximum parallel performance
        let mut data = Vec::with_capacity(file_size as usize);
        file.read_to_end(&mut data)?;

        if data.len() < 4 {
            return Ok(false);
        }

        let mut cache: BamCache;

        if &data[0..4] == CACHE_MAGIC_V1 {
            // New multi-block format
            let mut cursor = std::io::Cursor::new(&data[4..]);
            
            // 1. Read header compression block info (fast metadata)
            let compressed_headers: Option<crate::compression::CompressedHeaders> = bincode::deserialize_from(&mut cursor)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
            
            let zstd_pos = 4 + cursor.position() as usize;
            let zstd_data = &data[zstd_pos..];

            // If file is large enough, parallelize the data and header decompression tasks
            if file_size > 1024 * 1024 {
                let (header_res, cache_res) = rayon::join(
                    || {
                        if let Some(ref h) = compressed_headers {
                            let t_h = Instant::now();
                            let names = crate::compression::decompress_header_block(h);
                            Some((names, t_h.elapsed()))
                        } else {
                            None
                        }
                    },
                    || {
                        let decoder = zstd::stream::Decoder::new(zstd_data)?;
                        let res: BamCache = bincode::deserialize_from(decoder)
                            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
                        Ok::<BamCache, std::io::Error>(res)
                    }
                );

                cache = cache_res?;
                if let Some((names, elapsed)) = header_res {
                    cache.names = names;
                    println!("[BENCHMARK] Rust Manager: (Parallel) Header decompression took {:?}", elapsed);
                }
            } else {
                // Sequential for small files
                let decoder = zstd::stream::Decoder::new(zstd_data)?;
                cache = bincode::deserialize_from(decoder)
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
                
                if let Some(ref h) = compressed_headers {
                    let t_h = Instant::now();
                    cache.names = crate::compression::decompress_header_block(h);
                    println!("[BENCHMARK] Rust Manager: Header decompression took {:?}", t_h.elapsed());
                }
            }
        } else {
            // Legacy format or non-v1 file
            let decoder = zstd::stream::Decoder::new(&data[..])?;
            cache = bincode::deserialize_from(decoder)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
            
            // Fallback: Check if there's a compressed header block inside the cache object
            if let Some(ref compressed) = cache.compressed_headers {
                cache.names = crate::compression::decompress_header_block(compressed);
                cache.compressed_headers = None;
            }
        }

        let t_deser = Instant::now();
        println!("[BENCHMARK] Rust Manager: multi-block load for {} took {:?}", bam_path, t_deser.duration_since(t_start));

        // Re-calculate size_bytes (estimate)
        let t_size = Instant::now();
        let mut bytes = cache.reads.capacity() * std::mem::size_of::<CompactRead>();
        bytes += cache.segments.capacity() * std::mem::size_of::<CompactSegment>();
        bytes += cache.segment_tags.capacity();
        bytes += cache.names.capacity() * std::mem::size_of::<String>();
        for name in &cache.names {
            bytes += name.capacity();
        }
        bytes += cache.ref_names.capacity() * std::mem::size_of::<String>();
        for name in &cache.ref_names {
            bytes += name.capacity();
        }
        bytes += cache.ref_lengths.capacity() * std::mem::size_of::<usize>();
        
        cache.size_bytes = bytes;
        println!("[BENCHMARK] Rust Manager: Memory estimation took {:?}", t_size.elapsed());

        let normalized = normalize_path(bam_path);
        {
            let mut caches = self.caches.write().map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
            caches.insert(normalized, cache);
        }
        
        println!("[BENCHMARK] Rust Manager: Total disk-to-RAM load took {:?}", t_start.elapsed());
        self.benchmark.record_stage(&format!("Load Cache Disk: {}", bam_path), t_start.elapsed());
        self.benchmark.observe_peak();
        
        self.print_cache_layout();
        Ok(true)
    }


    pub fn get_or_load_compact_cache(&self, bam_path: &str) -> Option<(Vec<CompactRead>, Vec<CompactSegment>, Vec<u8>, Vec<String>, Vec<String>, Vec<usize>, CacheStatus)> {
        let normalized = normalize_path(bam_path);
        
        // Try RAM first
        {
            let caches = self.caches.read().ok()?;
            if let Some(cache) = caches.get(&normalized) {
                return Some((
                    cache.reads.clone(),
                    cache.segments.clone(),
                    cache.segment_tags.clone(),
                    cache.names.clone(),
                    cache.ref_names.clone(),
                    cache.ref_lengths.clone(),
                    cache.status
                ));
            }
        }

        // Try load from disk
        if self.load_from_disk(bam_path).unwrap_or(false) {
            let caches = self.caches.read().ok()?;
            if let Some(cache) = caches.get(&normalized) {
                return Some((
                    cache.reads.clone(),
                    cache.segments.clone(),
                    cache.segment_tags.clone(),
                    cache.names.clone(),
                    cache.ref_names.clone(),
                    cache.ref_lengths.clone(),
                    cache.status
                ));
            }
        }

        None
    }

    pub fn get_reads(&self, bam_path: &str) -> Option<(Vec<ReadInfo>, CacheStatus)> {
        // This is a bridge for the old API if needed, but we should use get_filtered_reads
        // to avoid reconstructing all reads in memory.
        let normalized = normalize_path(bam_path);
        let caches = self.caches.read().ok()?;
        let cache = caches.get(&normalized)?;
        let names = cache.materialized_names();
        
        if cache.status != CacheStatus::Complete {
            return Some((Vec::new(), cache.status));
        }

        // Warning: This reconstructs EVERYTHING. Only use if strictly necessary.
        let mut results = Vec::with_capacity(cache.reads.len());
        for read in &cache.reads {
            let mut segs = Vec::with_capacity(read.seg_count as usize);
            let is_mate = (read.flags & 0x02) != 0;
            for i in 0..read.seg_count as usize {
                let s = &cache.segments[read.seg_offset as usize + i];
                let tag = cache.segment_tags[read.seg_offset as usize + i];
                segs.push(ReadSegment { 
                    start: s.start as usize, 
                    end: s.end as usize, 
                    is_mate,
                    mismatches: (tag >> 1) & 0x0F,
                    insertions: (tag >> 5) & 0x07,
                    is_followed_by_deletion: (tag & 0x01) != 0,
                });
            }

            results.push(ReadInfo {
                name: names[read.name_id as usize].clone(),
                reference: cache.ref_names[read.ref_id as usize].clone(),
                mapping_quality: read.mq,
                strand: decode_strand_bits(read.flags),
                start: read.start as usize,
                end: read.end as usize,
                segments: segs,
            });
        }
        
        Some((results, cache.status))
    }

    pub fn get_filtered_reads(
        &self, 
        bam_path: &str, 
        ref_name: &str, 
        start_target: usize, 
        end_target: usize, 
        min_mq: u8
    ) -> Option<(Vec<ReadInfo>, CacheStatus)> {
        let normalized = normalize_path(bam_path);
        let caches = self.caches.read().ok()?;
        let cache = caches.get(&normalized)?;
        let names = cache.materialized_names();
        
        if cache.status != CacheStatus::Complete {
            return Some((Vec::new(), cache.status));
        }

        let target_ref_id = cache.ref_names.iter().position(|r| r == ref_name)? as u16;
        
        // Pass 1: Identify Name IDs that pass the junction or are parts of a bridge
        let mut passing_name_ids = std::collections::HashSet::new();
        let mut anchored_start = std::collections::HashSet::new();
        let mut anchored_end = std::collections::HashSet::new();

        for read in &cache.reads {
            if read.ref_id != target_ref_id || read.mq < min_mq { continue; }
            
            let mut has_junction = false;
            let mut has_start = false;
            let mut has_end = false;

            for i in 0..read.seg_count as usize {
                let seg = &cache.segments[read.seg_offset as usize + i];
                if i < (read.seg_count as usize).saturating_sub(1) {
                    let next_seg = &cache.segments[read.seg_offset as usize + i + 1];
                    let tag = cache.segment_tags[read.seg_offset as usize + i];
                    let is_deletion = (tag & 0x01) != 0;
                    if !is_deletion && seg.end as usize == start_target && next_seg.start as usize == end_target {
                        has_junction = true;
                    }
                }
                if seg.end as usize == start_target { has_start = true; }
                if seg.start as usize == end_target { has_end = true; }
            }

            if has_junction {
                passing_name_ids.insert(read.name_id);
            } else {
                if has_start { anchored_start.insert(read.name_id); }
                if has_end { anchored_end.insert(read.name_id); }
            }
        }

        // Bridge logic: if same name ID has one record at start and another at end
        for name_id in anchored_start {
            if anchored_end.contains(&name_id) {
                passing_name_ids.insert(name_id);
            }
        }

        // Pass 2: Reconstruct ReadInfo only for passing names
        let mut results = Vec::new();
        for read in &cache.reads {
            if passing_name_ids.contains(&read.name_id) {
                let mut segs = Vec::with_capacity(read.seg_count as usize);
                let is_mate = (read.flags & 0x02) != 0;
                for i in 0..read.seg_count as usize {
                    let s = &cache.segments[read.seg_offset as usize + i];
                    let tag = cache.segment_tags[read.seg_offset as usize + i];
                    segs.push(ReadSegment { 
                        start: s.start as usize, 
                        end: s.end as usize, 
                        is_mate,
                        mismatches: (tag >> 1) & 0x0F,
                        insertions: (tag >> 5) & 0x07,
                        is_followed_by_deletion: (tag & 0x01) != 0,
                    });
                }

                results.push(ReadInfo {
                    name: names[read.name_id as usize].clone(),
                    reference: ref_name.to_string(),
                    mapping_quality: read.mq,
                    strand: decode_strand_bits(read.flags),
                    start: read.start as usize,
                    end: read.end as usize,
                    segments: segs,
                });
            }
        }

        Some((results, cache.status))
    }

    pub fn get_filtered_reads_batch(
        &self,
        bam_path: &str,
        ref_name: &str,
        junctions: &[(usize, usize)],
        min_mq: u8,
    ) -> Option<(std::collections::HashMap<(usize, usize), Vec<ReadInfo>>, CacheStatus)> {
        let normalized = normalize_path(bam_path);
        let caches = self.caches.read().ok()?;
        let cache = caches.get(&normalized)?;
        let names = cache.materialized_names();

        if cache.status != CacheStatus::Complete {
            return Some((std::collections::HashMap::new(), cache.status));
        }

        let target_ref_id = cache.ref_names.iter().position(|r| r == ref_name)? as u16;

        // Map junctions to sets of passing name IDs
        let mut junc_to_passing_names = std::collections::HashMap::new();
        let mut junc_to_anchored_start = std::collections::HashMap::new();
        let mut junc_to_anchored_end = std::collections::HashMap::new();

        for &(js, je) in junctions {
            junc_to_passing_names.insert((js, je), std::collections::HashSet::new());
            junc_to_anchored_start.insert((js, je), std::collections::HashSet::new());
            junc_to_anchored_end.insert((js, je), std::collections::HashSet::new());
        }

        // Global map to quickly find all junctions associated with a name ID
        let mut name_id_to_juncs = std::collections::HashMap::new();

        // Single pass over all reads
        for read in &cache.reads {
            if read.ref_id != target_ref_id || read.mq < min_mq {
                continue;
            }

            for &(js, je) in junctions {
                let mut has_junction = false;
                let mut has_start = false;
                let mut has_end = false;

                for i in 0..read.seg_count as usize {
                    let seg = &cache.segments[read.seg_offset as usize + i];
                    if i < (read.seg_count as usize).saturating_sub(1) {
                        let next_seg = &cache.segments[read.seg_offset as usize + i + 1];
                        let tag = cache.segment_tags[read.seg_offset as usize + i];
                        let is_deletion = (tag & 0x01) != 0;
                        if !is_deletion && seg.end as usize == js && next_seg.start as usize == je {
                            has_junction = true;
                        }
                    }
                    if seg.end as usize == js {
                        has_start = true;
                    }
                    if seg.start as usize == je {
                        has_end = true;
                    }
                }

                if has_junction {
                    junc_to_passing_names.get_mut(&(js, je)).unwrap().insert(read.name_id);
                    name_id_to_juncs.entry(read.name_id).or_insert_with(std::collections::HashSet::new).insert((js, je));
                } else {
                    if has_start {
                        junc_to_anchored_start.get_mut(&(js, je)).unwrap().insert(read.name_id);
                    }
                    if has_end {
                        junc_to_anchored_end.get_mut(&(js, je)).unwrap().insert(read.name_id);
                    }
                }
            }
        }

        // Bridge logic for each junction
        for (&(js, je), anchored_start) in &junc_to_anchored_start {
            let anchored_end = junc_to_anchored_end.get(&(js, je)).unwrap();
            let passing_names = junc_to_passing_names.get_mut(&(js, je)).unwrap();
            for &name_id in anchored_start {
                if anchored_end.contains(&name_id) {
                    passing_names.insert(name_id);
                    name_id_to_juncs.entry(name_id).or_insert_with(std::collections::HashSet::new).insert((js, je));
                }
            }
        }

        // Final reconstruction: Pass 2
        let mut batch_results = std::collections::HashMap::new();
        for &j in junctions {
            batch_results.insert(j, Vec::new());
        }

        for read in &cache.reads {
            if let Some(matching_juncs) = name_id_to_juncs.get(&read.name_id) {
                let mut segs = Vec::with_capacity(read.seg_count as usize);
                let is_mate = (read.flags & 0x02) != 0;
                for i in 0..read.seg_count as usize {
                    let s = &cache.segments[read.seg_offset as usize + i];
                    let tag = cache.segment_tags[read.seg_offset as usize + i];
                    segs.push(ReadSegment {
                        start: s.start as usize,
                        end: s.end as usize,
                        is_mate,
                        mismatches: (tag >> 1) & 0x0F,
                        insertions: (tag >> 5) & 0x07,
                        is_followed_by_deletion: (tag & 0x01) != 0,
                    });
                }

                let info = ReadInfo {
                    name: names[read.name_id as usize].clone(),
                    reference: ref_name.to_string(),
                    mapping_quality: read.mq,
                    strand: decode_strand_bits(read.flags),
                    start: read.start as usize,
                    end: read.end as usize,
                    segments: segs,
                };

                for &j in matching_juncs {
                    batch_results.get_mut(&j).unwrap().push(info.clone());
                }
            }
        }

        Some((batch_results, cache.status))
    }

    pub fn start_caching(&self, bam_path: &str) {
        let normalized = normalize_path(bam_path);
        if let Ok(mut caches) = self.caches.write() {
            caches.insert(normalized, BamCache {
                reads: Vec::new(),
                segments: Vec::new(),
                segment_tags: Vec::new(),
                names: Vec::new(),
                delta_names: None,
                ref_names: Vec::new(),
                ref_lengths: Vec::new(),
                status: CacheStatus::Caching,
                compressed_headers: None,
                size_bytes: 0,
                zstd_payload_uncompressed: 0,
                zstd_payload_compressed: 0,
            });
        }
    }

    pub fn commit_caching_compact(
        &self, 
        bam_path: &str, 
        raw_reads: Vec<ReadInfo>, 
        ref_names: Vec<String>,
        ref_lengths: Vec<usize>,
        status: CacheStatus,
    ) -> std::io::Result<()> {
        let normalized = normalize_path(bam_path);
        
        if status == CacheStatus::ExceededLimit {
            return Err(std::io::Error::new(std::io::ErrorKind::Other, "Cache status is ExceededLimit"));
        }

        let provider = get_read_provider();
        let use_delta = *provider.use_delta_incremental.read().unwrap_or_else(|e| e.into_inner());
        let use_substitutes = *provider.use_substitutes.read().unwrap_or_else(|e| e.into_inner());

        let mut name_to_id = HashMap::with_capacity(if use_delta { 0 } else { raw_reads.len() / 2 });
        let mut name_hash_to_id = HashMap::with_capacity(if use_delta { raw_reads.len() / 2 } else { 0 });
        
        let mut names = Vec::new();
        let mut delta_names = if use_delta { Some(DeltaNameStore::new(use_substitutes)) } else { None };
        
        let mut segments = Vec::new();
        let mut segment_tags = Vec::new();
        let mut compact_reads = Vec::with_capacity(raw_reads.len());

        let ref_to_id: HashMap<String, u16> = ref_names.iter().enumerate()
            .map(|(i, n)| (n.clone(), i as u16)).collect();

        for read in raw_reads {
            let name_id = if let Some(ref mut delta) = delta_names {
                let hash = hash_name_128(&read.name);
                *name_hash_to_id.entry(hash).or_insert_with(|| {
                    let id = delta.num_names as u32;
                    delta.add_name(&read.name);
                    id
                })
            } else {
                *name_to_id.entry(read.name.clone()).or_insert_with(|| {
                    let id = names.len() as u32;
                    names.push(read.name);
                    id
                })
            };

            let ref_id = *ref_to_id.get(&read.reference).unwrap_or(&0);
            
            // Flags: bit 0 = paired, bit 1 = is_mate
            let mut flags = 0u8;
            flags |= 0x01; // Assume aligned records are paired in this context
            if read.segments.first().map(|s| s.is_mate).unwrap_or(false) {
                flags |= 0x02;
            }
            flags |= encode_strand_bits(read.strand);

            let seg_offset = segments.len() as u32;
            let seg_count = read.segments.len() as u8;
            for s in read.segments {
                segments.push(CompactSegment { start: s.start as u32, end: s.end as u32 });
                
                let mut tag = 0u8;
                if s.is_followed_by_deletion { tag |= 0x01; }
                tag |= (s.mismatches & 0x0F) << 1;
                tag |= (s.insertions & 0x07) << 5;
                segment_tags.push(tag);
            }

            compact_reads.push(CompactRead {
                name_id,
                ref_id,
                mq: read.mapping_quality,
                flags,
                start: read.start as u32,
                end: read.end as u32,
                seg_offset,
                seg_count,
            });
        }

        // Calculate size carefully including capacities
        let mem_reads = compact_reads.capacity() * std::mem::size_of::<CompactRead>();
        let mem_segs = segments.capacity() * std::mem::size_of::<CompactSegment>();
        let mem_tags = segment_tags.capacity();
        let mem_names: usize = names.capacity() * std::mem::size_of::<String>() + 
                               names.iter().map(|s| s.capacity()).sum::<usize>();
        let mem_ref_names: usize = ref_names.capacity() * std::mem::size_of::<String>() + 
                                   ref_names.iter().map(|s| s.capacity()).sum::<usize>();
        let mem_ref_lengths = ref_lengths.capacity() * std::mem::size_of::<usize>();
        
        let size_bytes = mem_reads + mem_segs + mem_tags + mem_names + mem_ref_names + mem_ref_lengths;
        let mem_mb = size_bytes as f64 / 1024.0 / 1024.0;

        // Record BAM stats
        let mapping_mem_id = name_to_id.capacity() * (std::mem::size_of::<String>() + std::mem::size_of::<u32>()) +
                             name_to_id.keys().map(|s| s.capacity()).sum::<usize>();
        let mapping_mem_hash = name_hash_to_id.capacity() * (std::mem::size_of::<u128>() + std::mem::size_of::<u32>());
        let mapping_mem_ref = ref_to_id.capacity() * (std::mem::size_of::<String>() + std::mem::size_of::<u16>()) +
                              ref_to_id.keys().map(|s| s.capacity()).sum::<usize>();
        let mapping_mb = (mapping_mem_id + mapping_mem_hash + mapping_mem_ref) as f64 / 1024.0 / 1024.0;
        let os_mem = crate::benchmarking::get_os_memory();

        self.benchmark.record_bam_stats(crate::benchmarking::BamStats {
            context: normalized.clone(),
            read_count: compact_reads.len(),
            segment_count: segments.len(),
            data_mb: mem_mb,
            mapping_mb,
            os_mb: os_mem as f64 / 1024.0 / 1024.0,
            header_uncompressed_bytes: names.iter().map(|s| s.len() as u64).sum(),
            header_compressed_bytes: 0, // Not compressed yet
            payload_uncompressed_bytes: size_bytes as u64 - (mem_names + mem_ref_names + mem_ref_lengths) as u64,
            payload_compressed_bytes: size_bytes as u64 - (mem_names + mem_ref_names + mem_ref_lengths) as u64,
        });

        let current_usage = self.get_total_memory_usage_mb();
        let limit = *self.max_memory_mb.read().unwrap_or_else(|e| e.into_inner());

        if current_usage + mem_mb > limit {
            let msg = format!("Rust Cache: Memory limit exceeded ({:.1} MB usage + {:.1} MB new > {:.1} MB limit) for {}.", 
                     current_usage, mem_mb, limit, normalized);
            println!("{}", msg);
            return Err(std::io::Error::new(std::io::ErrorKind::Other, msg));
        }

        if let Ok(mut caches) = self.caches.write() {
            let mode_str = if use_delta { "Delta-Incremental" } else { "Standard" };
            println!("Rust: Collection complete ({} mode) for {}: {} reads, {} segments, {} unique names (~{:.2} MB)", 
                     mode_str, normalized, compact_reads.len(), segments.len(), 
                     if let Some(ref d) = delta_names { d.num_names } else { names.len() }, 
                     mem_mb);
            
            let mut cache = BamCache {
                reads: compact_reads,
                segments,
                segment_tags,
                names,
                delta_names,
                ref_names,
                ref_lengths,
                status,
                compressed_headers: None,
                size_bytes,
                zstd_payload_uncompressed: 0,
                zstd_payload_compressed: 0,
            };

            if self.compression_mode.read().unwrap().is_some() && status == CacheStatus::Complete {
                cache.compress_headers();
            }

            caches.insert(normalized.clone(), cache);
        }
        
        // Persist to disk if complete
        if status == CacheStatus::Complete {
            let _ = self.persist_to_disk(bam_path);
        }

        self.print_cache_layout();
        Ok(())
    }

    pub fn commit_builder(&self, bam_path: &str, builder: IncrementalCacheBuilder, status: CacheStatus) -> std::io::Result<()> {
        let normalized = normalize_path(bam_path);
        
        let mapping_mb = builder.get_mapping_memory_estimate() as f64 / 1024.0 / 1024.0;
        let read_count = builder.cache.reads.len();
        let segment_count = builder.cache.segments.len();
        let os_mem = crate::benchmarking::get_os_memory();
        
        let mut cache = builder.finalize(status);
        let data_mb = cache.size_bytes as f64 / 1024.0 / 1024.0;

        self.benchmark.record_bam_stats(crate::benchmarking::BamStats {
            context: normalized.clone(),
            read_count,
            segment_count,
            data_mb,
            mapping_mb,
            os_mb: os_mem as f64 / 1024.0 / 1024.0,
            header_uncompressed_bytes: cache.get_uncompressed_header_bytes(),
            header_compressed_bytes: cache.get_compressed_header_bytes(),
            payload_uncompressed_bytes: cache.get_payload_bytes(),
            payload_compressed_bytes: cache.get_payload_bytes(),
        });

        let mem_mb = data_mb;
        let current_usage = self.get_total_memory_usage_mb();
        let limit = *self.max_memory_mb.read().unwrap_or_else(|e| e.into_inner());

        if current_usage + mem_mb > limit {
            let msg = format!("Rust Cache: Memory limit exceeded ({:.1} MB usage + {:.1} MB new > {:.1} MB limit) for {}.", 
                     current_usage, mem_mb, limit, normalized);
            println!("{}", msg);
            return Err(std::io::Error::new(std::io::ErrorKind::Other, msg));
        }

        if let Ok(mut caches) = self.caches.write() {
            let mode_str = if cache.delta_names.is_some() { "Delta-Incremental" } else { "Standard" };
            println!("Rust: Collection complete ({} mode) for {}: {} reads, {} segments (~{:.2} MB)", 
                     mode_str, normalized, cache.reads.len(), cache.segments.len(), mem_mb);

            if self.compression_mode.read().unwrap().is_some() && status == CacheStatus::Complete {
                cache.compress_headers();
            }

            caches.insert(normalized.clone(), cache);
        }

        if status == CacheStatus::Complete {
            let _ = self.persist_to_disk(bam_path);
        }

        self.print_cache_layout();
        Ok(())
    }

    pub fn commit_caching(&self, _bam_path: &str, _reads: Vec<ReadInfo>, _status: CacheStatus) {
        // Obsolete, replaced by commit_caching_compact
    }

    pub fn clear_cache(&self) {
        if let Ok(mut caches) = self.caches.write() {
            caches.clear();
        }
        self.print_cache_layout();
    }

    pub fn estimate_memory_mb(&self, raw_reads: &[ReadInfo]) -> f64 {
        let mem_reads = raw_reads.len() * std::mem::size_of::<CompactRead>();
        let mem_segs: usize = raw_reads.iter().map(|r| r.segments.len() * std::mem::size_of::<CompactSegment>()).sum();
        // Use 48 bytes per string for name overhead (Vec pointer + String pointer/len/cap + allocation overhead)
        let mem_names: usize = raw_reads.iter().map(|r| r.name.len() + 48).sum();
        ((mem_reads + mem_segs + mem_names) as f64 / 1024.0 / 1024.0) * 1.1 // 10% safety margin
    }
}

#[cfg(test)]
mod manager_tests {
    use super::*;

    #[test]
    fn test_delta_incremental_collection() {
        let provider = get_read_provider();
        provider.set_delta_incremental(true);
        provider.set_compress_headers(Some(CompressionMode::Huffman), true);
        
        let ref_names = vec!["chr1".to_string()];
        let ref_lengths = vec![1000];
        let mut builder = IncrementalCacheBuilder::new("test.bam", ref_names.clone(), ref_lengths.clone());
        
        let read1 = ReadInfo {
            name: "INSTRUMENT:1:FLOWCELL:1:1001".to_string(),
            reference: "chr1".to_string(),
            mapping_quality: 60,
            strand: 0,
            start: 100,
            end: 200,
            segments: vec![ReadSegment { start: 100, end: 200, is_mate: false, mismatches: 0, insertions: 0, is_followed_by_deletion: false }],
        };
        
        let read2 = ReadInfo {
            name: "INSTRUMENT:1:FLOWCELL:1:1002".to_string(),
            reference: "chr1".to_string(),
            mapping_quality: 60,
            strand: 0,
            start: 150,
            end: 250,
            segments: vec![ReadSegment { start: 150, end: 250, is_mate: false, mismatches: 0, insertions: 0, is_followed_by_deletion: false }],
        };
        
        builder.add_read(read1.clone());
        builder.add_read(read2.clone());
        
        let cache = builder.finalize(CacheStatus::Complete);
        assert!(cache.delta_names.is_some());
        assert_eq!(cache.delta_names.as_ref().unwrap().num_names, 2);
        
        let names = cache.delta_names.as_ref().unwrap().to_names();
        assert_eq!(names.len(), 2);
        assert_eq!(names[0], read1.name);
        assert_eq!(names[1], read2.name);
        
        // Test transition to compressed headers
        let mut cache_final = cache;
        cache_final.compress_headers();
        assert!(cache_final.compressed_headers.is_some());
        assert!(cache_final.delta_names.is_none());
        assert!(cache_final.names.is_empty());
        
        // Verify we can still get the names back via decompression
        let decompressed = crate::compression::decompress_header_block(cache_final.compressed_headers.as_ref().unwrap());
        assert_eq!(decompressed.len(), 2);
        assert_eq!(decompressed[0], read1.name);
        assert_eq!(decompressed[1], read2.name);
    }
}

