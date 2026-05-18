// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
//! Pipeline entrypoint: parse config, compute introns, collect reads, write BAM+BAI

mod extract;
mod gtf_query;
mod collect;
mod config;
mod header_map;
mod intron;
mod io;

use noodles::bam;
use noodles::core::Region;
use std::collections::{HashSet, HashMap};
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::fs;
use serde_json::json;
use clap::Parser;

use collect::{collect_records_from_bams, collect_records_streamed};
use config::{AnyResult, parse_config};
use lnc_seeker_lib::{AssemblyReport, ProgressData};
use lnc_seeker_lib::reads_manager::{get_read_provider, BamCache};
use lnc_seeker_lib::compression::CompressionMode;
use lnc_seeker_lib::analysis::{populate_bam_cache, populate_cache_from_records};
use lnc_seeker_lib::genome::GenomeProvider;
use header_map::build_name_to_index;
use intron::get_gene_regions;
use io::{sort_records_for_indexing, write_bam_and_bai};

#[derive(Parser, Debug)]
#[command(author, version, about = "lnc-seeker-collect: collect reads from BAM files based on regions and cohorts")]
struct Args {
    /// Path to the main configuration file
    #[arg(short, long, default_value = "config.cfg")]
    config: String,

    /// Override the GTF annotation file path
    #[arg(short, long)]
    gtf: Option<String>,

    /// Override the assembly report (chromosome mapping) file path
    #[arg(short, long)]
    assembly_report: Option<String>,

    /// Override the genome FASTA file path
    #[arg(long)]
    genome: Option<String>,

    /// Override the output directory for BAM files and dictionary
    #[arg(short, long)]
    output_dir: Option<String>,

    /// Override the path prefix for BAM paths in dictionary.json
    #[arg(short = 'P', long)]
    path_prefix: Option<String>,

    /// Additional configuration files to include and parse
    #[arg(short, long = "include")]
    extra_includes: Vec<String>,

    /// Override gene regions (replaces all defined in config). 
    /// Format: "gene1,gene2[:offset]" (repeatable)
    #[arg(short = 'G', long = "gene-region")]
    gene_regions: Vec<String>,

    /// Override extra regions (replaces all defined in config). 
    /// Format: "chr1:100-200" (repeatable)
    #[arg(short = 'R', long = "region")]
    extra_regions: Vec<String>,

    /// Keep secondary alignments (they are ignored by default)
    #[arg(long, default_value_t = false)]
    keep_secondary: bool,

    /// Force endpoint annotation tag (EB/EE) enablement
    #[arg(long, default_value_t = false)]
    annotate_endpoint: bool,

    /// Force recreation of BAM files even if they already exist
    #[arg(short, long, default_value_t = false)]
    force: bool,

    /// Use incremental streaming to avoid memory spikes (experimental, default)
    #[arg(long)]
    stream: bool,

    /// Disable incremental streaming (force traditional batch mode)
    #[arg(long)]
    no_stream: bool,

    /// Verify streaming results against batch results (very memory intensive)
    #[arg(long, default_value_t = false)]
    verify: bool,

    /// Use experimental differential Huffman header compression (saves disk space)
    #[arg(long, default_value_t = false)]
    compress_headers: bool,

    /// Use incremental name storage to reduce memory peaks during collection
    #[arg(long)]
    use_delta_collection: bool,

    /// Algorithm for header compression: "huffman", "zstd", or "delta" (default: "huffman")
    #[arg(long)]
    header_compression_algo: Option<String>,
}

fn main() -> AnyResult<()> {
    let args = Args::parse();

    println!("lnc-seeker-collect: Loading configuration from {}", args.config);

    let mut config = parse_config(&args.config, &args.extra_includes)?;

    // Apply CLI overrides
    if let Some(gtf) = args.gtf {
        config.gtf_path = gtf;
    }
    if let Some(ar) = args.assembly_report {
        config.assembly_report = ar;
    }
    if let Some(genome) = args.genome {
        config.genome_path = Some(genome);
    }
    if let Some(od) = args.output_dir {
        config.output_dir = Some(od);
    }
    if let Some(pp) = args.path_prefix {
        config.path_prefix = Some(pp);
    }
    if args.keep_secondary {
        config.skip_secondary = false;
    }

    if config.genome_path.is_none() {
        eprintln!("Warning: No genome path provided. Compact storage caches will be generated without mismatch information.");
    } else if let Some(ref path) = config.genome_path {
        if !Path::new(path).exists() {
            eprintln!("Warning: Genome file '{}' not found. Compact storage caches will be generated without mismatch information.", path);
        }
    }

    if args.annotate_endpoint {
        config.annotate_endpoint_tag = true;
    }

    if args.compress_headers {
        config.compress_headers = true;
    }

    if args.use_delta_collection {
        config.use_delta_incremental_collection = true;
    }

    if let Some(algo) = args.header_compression_algo {
        config.compression_algorithm = algo.to_ascii_lowercase();
    }

    if args.stream {
        config.stream = true;
    }
    if args.no_stream {
        config.stream = false;
    }

    // Replace gene regions if specified on CLI
    if !args.gene_regions.is_empty() {
        config.gene_regions.clear();
        let mut seen_genes = HashSet::new();
        for val in &args.gene_regions {
            let (genes_part, offset) = if let Some((g, o)) = val.split_once(':') {
                let offset = o.trim().parse::<usize>()
                    .map_err(|e| format!("invalid CLI offset '{}': {}", o, e))?;
                (g, offset)
            } else {
                (val.as_str(), 0)
            };

            let mut line_genes = Vec::new();
            for gene in genes_part.split(',') {
                let gene = gene.trim().to_string();
                if gene.is_empty() { continue; }
                if !seen_genes.insert(gene.clone()) {
                    return Err(format!("duplicate gene name in CLI arguments: {}", gene).into());
                }
                line_genes.push(gene);
            }
            if !line_genes.is_empty() {
                config.gene_regions.push((line_genes, offset));
            }
        }
    }

    // Replace extra regions if specified on CLI
    if !args.extra_regions.is_empty() {
        config.extra_regions = args.extra_regions;
    }

    let assembly_report =
        AssemblyReport::parse_from_file(&config.assembly_report)?;

    // Collect all genes across all gene_region lines to fetch their data in one pass
    let all_genes: Vec<String> = config.gene_regions.iter().flat_map(|(genes, _)| genes.clone()).collect();
    let gene_data_map = get_gene_regions(&config.gtf_path, &assembly_report, &all_genes)?;

    // 1. Compute all extraction regions indicated by the various parameters in the config-file.
    let mut regions = Vec::new();
    let mut dict_metadata = Vec::new();

    // Extra regions from config
    for r_str in &config.extra_regions {
        let region: Region = r_str
            .parse()
            .map_err(|e| format!("invalid region '{}': {}", r_str, e))?;
        regions.push(region.clone());
        dict_metadata.push((r_str.clone(), 0, region.to_string()));
    }

    // Gene regions (one per line in config)
    for (line_genes, offset) in &config.gene_regions {
        let gene_region = extract::get_overall_gene_region(line_genes.clone(), &gene_data_map, *offset)?;
        regions.push(gene_region.clone());
        dict_metadata.push((line_genes.join(","), *offset, gene_region.to_string()));
    }

    // 2. Use the regions computed in 1. to collect the reads via using an instance of RegionSet.
    let region_set = lnc_seeker_lib::regions::RegionSet::new(regions, &assembly_report);
    region_set.pretty_print();

    // Generate unique tags *per region* to avoid filename conflicts.
    // Each region entry in the final dictionary will point to its own set of BAMs.
    let mut region_tags = Vec::new();
    for (key, offset, region_str) in &dict_metadata {
        let mut h = DefaultHasher::new();
        key.hash(&mut h);
        offset.hash(&mut h);
        region_str.hash(&mut h);
        let hash_str = format!("{:08x}", h.finish());
        let prefix = key.split(',').next().unwrap_or("multi")
            .chars().filter(|c| c.is_alphanumeric()).collect::<String>();
        region_tags.push(format!("{}_{}", prefix, hash_str));
    }

    // Map: region_key -> Map: cohort_name -> cohort_info
    let mut region_to_cohort_map: HashMap<String, serde_json::Map<String, serde_json::Value>> = HashMap::new();

    let _out_dir = if let Some(ref dir) = config.output_dir {
        let out_dir = Path::new(dir);
        if !out_dir.exists() {
            fs::create_dir_all(out_dir).map_err(|e| format!("Failed to create output directory {}: {}", dir, e))?;
        }
        Some(out_dir)
    } else {
        None
    };

    let read_provider = get_read_provider();
    // Set a very high memory limit for the collector, as it's a batch process 
    // and needs to be able to handle large regions like MALAT1.
    read_provider.set_max_memory_mb(1_048_576.0); // 1 TB
    let compression_mode = if config.compress_headers {
        match config.compression_algorithm.as_str() {
            "zstd" => Some(CompressionMode::Zstd),
            "delta" | "none" => Some(CompressionMode::None),
            _ => Some(CompressionMode::Huffman),
        }
    } else {
        None
    };
    read_provider.set_compress_headers(compression_mode, config.compression_use_substitutes);
    read_provider.set_delta_incremental(config.use_delta_incremental_collection);
    
    let progress = Arc::new(ProgressData {
        stage: std::sync::atomic::AtomicU32::new(0),
        current: std::sync::atomic::AtomicU32::new(0),
        total: std::sync::atomic::AtomicU32::new(0),
    });

    let mut genome_provider = GenomeProvider::new(&config.genome_path)?;

    for cohort in &config.cohorts {
        println!("\n=== Processing Cohort: {} ===", cohort.name);
        println!("Input BAMs: {} file(s)", cohort.bam_paths.len());

        read_provider.benchmark.set_current_cohort(Some(cohort.name.clone()));
        let t_cohort_begin = std::time::Instant::now();

        // Read header from the first input BAM to preserve reference sequences.
        let first_bam = cohort.bam_paths.first().expect("no bam paths");
        if !std::path::Path::new(first_bam).exists() {
            return Err(format!("First BAM of cohort '{}' not found: {}", cohort.name, first_bam).into());
        }

        let mut out_ireader = bam::io::indexed_reader::Builder::default()
            .build_from_path(first_bam)
            .map_err(|e| format!("{}: {}", first_bam, e))?;
        let out_header = out_ireader
            .read_header()
            .map_err(|e| format!("{}: {}", first_bam, e))?;
        let out_name_to_index = build_name_to_index(&out_header);

        // 1. Resolve all output paths and check for existence
        let mut output_paths = Vec::new();
        let mut need_collection_indices = Vec::new();
        let mut skipped_count = 0;

        for (idx, _) in dict_metadata.iter().enumerate() {
            let region_tag = &region_tags[idx];
            let output_bam_path = if let Some(ref manual_path) = cohort.output_bam {
                let p = Path::new(manual_path);
                let stem = p.file_stem().map(|s| s.to_string_lossy()).unwrap_or_else(|| "output".into());
                let ext = p.extension().map(|e| e.to_string_lossy()).unwrap_or_else(|| "bam".into());
                p.with_file_name(format!("{}_{}.{}", stem, region_tag, ext)).to_string_lossy().to_string()
            } else {
                let dir = config.output_dir.as_ref().expect("output_dir or output_bam must be set");
                let mut path = PathBuf::from(dir);
                path.push(format!("{}_{}.bam", cohort.name, region_tag));
                path.to_string_lossy().to_string()
            };
            
            let cache_path = {
                let mut p = PathBuf::from(&output_bam_path);
                p.set_extension("lnc_cache.bin");
                p
            };

            let exists = if config.write_bam {
                Path::new(&output_bam_path).exists()
            } else {
                cache_path.exists()
            };
            
            if !args.force && exists {
                skipped_count += 1;
            } else {
                need_collection_indices.push(idx);
            }
            output_paths.push(output_bam_path);
        }

        if skipped_count > 0 {
            let reason = if config.write_bam { "BAM files already exist" } else { "compact storage caches already exist" };
            println!("Note: Skipped collection for {} regions in cohort {} because {}.", skipped_count, cohort.name, reason);
        }

        let mut streamed_caches: Vec<Option<BamCache>> = vec![None; region_set.len()];
        if (config.stream || args.verify) && !need_collection_indices.is_empty() {
            println!("Starting STREAMED collection for cohort {} ({} regions)...", cohort.name, need_collection_indices.len());
            let mut subset_region_set = region_set.clone();
            subset_region_set.filter_by_indices(&need_collection_indices);

            let streamed_results = collect_records_streamed(
                &cohort.bam_paths,
                &subset_region_set,
                config.skip_secondary,
                &out_header,
                &assembly_report,
                &mut genome_provider,
                &cohort.name,
            )?;

            let mut results_iter = streamed_results.into_iter();
            for &orig_idx in &need_collection_indices {
                streamed_caches[orig_idx] = results_iter.next();
            }
        }

        let mut grouped_collected = if (!config.stream || args.verify) && !need_collection_indices.is_empty() {
            println!("Starting BATCH collection for cohort {} ({} regions)...", cohort.name, need_collection_indices.len());
            let mut subset_region_set = region_set.clone();
            subset_region_set.filter_by_indices(&need_collection_indices);
            
            let mut results = collect_records_from_bams(
                &cohort.bam_paths,
                &subset_region_set,
                config.annotate_endpoint_tag,
                config.skip_secondary,
                &out_name_to_index,
                &assembly_report,
                &cohort.name,
            )?;
            
            // Map results back to original indices
            let mut full_results = vec![Vec::new(); region_set.len()];
            for (i, &orig_idx) in need_collection_indices.iter().enumerate() {
                full_results[orig_idx] = std::mem::take(&mut results[i]);
            }
            full_results
        } else {
            vec![Vec::new(); region_set.len()]
        };

        // For each region, write newly collected BAMs and update map for all valid ones
        for (idx, (key, _offset, _region_str)) in dict_metadata.iter().enumerate() {
            read_provider.benchmark.observe_cohort_peak(&cohort.name);
            let output_bam_path = &output_paths[idx];
            let mut collected = std::mem::take(&mut grouped_collected[idx]);
            
            let cache_path = {
                let mut p = PathBuf::from(output_bam_path);
                p.set_extension("lnc_cache.bin");
                p
            };
            
            let exists_already = if !args.force {
                if config.write_bam {
                    Path::new(output_bam_path).exists()
                } else {
                    cache_path.exists()
                }
            } else {
                false
            };

            if !exists_already {
                let mut has_data = false;
                if args.verify || !config.stream {
                    if !collected.is_empty() {
                        has_data = true;
                        sort_records_for_indexing(&mut collected);

                        if config.write_bam {
                            println!("Region {}: Writing {} unique reads to {}...", key, collected.len(), output_bam_path);
                            write_bam_and_bai(output_bam_path, &out_header, &collected)?;
                        }
                    }
                } else {
                    // Stream mode (no verify)
                    if let Some(ref cache) = streamed_caches[idx] {
                        if !cache.reads.is_empty() {
                            has_data = true;
                        }
                    }
                }

                if !has_data {
                    println!("Region {}: No reads found.", key);
                    continue;
                }

                if args.verify {
                    let streamed = std::mem::take(&mut streamed_caches[idx]).expect("Missing streamed cache in verify mode");
                    
                    // Run batch compaction
                    println!("Region {}: Running batch compaction for verification...", key);
                    populate_cache_from_records(
                        output_bam_path,
                        &collected,
                        &out_header,
                        &mut genome_provider,
                        read_provider,
                    ).map_err(|e| format!("Batch population failed during verify: {}", e))?;
                    
                    let batch = {
                        let normalized = lnc_seeker_lib::utils::normalize_path(output_bam_path);
                        let caches = read_provider.caches.read().unwrap();
                        caches.get(&normalized).cloned().expect("Batch cache not found after population")
                    };

                    // Compare
                    if streamed.reads.len() != batch.reads.len() || streamed.names.len() != batch.names.len() || streamed.segments.len() != batch.segments.len() {
                        println!("Region {}: VERIFICATION FAILED!", key);
                        println!("  Streamed: {} reads, {} names, {} segments", streamed.reads.len(), streamed.names.len(), streamed.segments.len());
                        println!("  Batch:    {} reads, {} names, {} segments", batch.reads.len(), batch.names.len(), batch.segments.len());
                    } else {
                        println!("Region {}: VERIFICATION PASSED. {} reads matched.", key, streamed.reads.len());
                    }
                } else if config.stream {
                    if let Some(cache) = std::mem::take(&mut streamed_caches[idx]) {
                        println!("Region {}: Writing streamed compact storage cache directly to disk...", key);
                        
                        let normalized = lnc_seeker_lib::utils::normalize_path(output_bam_path);
                        if let Ok(mut caches) = read_provider.caches.write() {
                            caches.insert(normalized.clone(), cache);
                        }
                        read_provider.persist_to_disk(output_bam_path)
                            .map_err(|e| format!("Failed to persist streamed cache for {}: {}", output_bam_path, e))?;
                    }
                } else {
                    // Pure Batch mode
                    println!("Region {}: Generating compact storage cache directly from memory...", key);
                    populate_cache_from_records(
                        output_bam_path,
                        &collected,
                        &out_header,
                        &mut genome_provider,
                        read_provider,
                    ).map_err(|e| format!("Failed to generate direct cache for {}: {}", output_bam_path, e))?;
                }
            } else {
                // Backward compatibility: If it already exists, ensure compact storage (cache) exists
                if !cache_path.exists() {
                    println!("Region {}: Generating compact storage cache from existing BAM...", key);
                    populate_bam_cache(
                        output_bam_path,
                        &progress,
                        read_provider,
                        &config.genome_path,
                    ).map_err(|e| format!("Failed to generate cache for {}: {}", output_bam_path, e))?;
                } else {
                    // Populate metadata from existing cache file for inclusion in dictionary.json
                    read_provider.load_from_disk(output_bam_path)
                        .map_err(|e| format!("Failed to load existing cache for {}: {}", output_bam_path, e))?;
                }
            }
            
            // Extract metadata from cache before clearing
            let (num_reads, cache_size, total_bases) = {
                let normalized = lnc_seeker_lib::utils::normalize_path(output_bam_path);
                if let Ok(caches) = read_provider.caches.read() {
                    caches.get(&normalized).map(|c| (c.reads.len(), c.size_bytes, c.total_coverage_bases())).unwrap_or((0, 0, 0))
                } else {
                    (0, 0, 0)
                }
            };

            // Calculate region length for coverage normalization
            let region = &region_set.regions()[idx];
            let region_length = {
                let interval = region.interval();
                let start = interval.start().map(usize::from).unwrap_or(1);
                let end = interval.end().map(usize::from).unwrap_or_else(|| {
                    let chrom_name = region.name();
                    out_header.reference_sequences()
                        .get(chrom_name)
                        .map(|rs| rs.length().get())
                        .unwrap_or(1)
                });
                if end >= start { (end - start + 1) as f64 } else { 1.0 }
            };

            let avg_coverage_per_sample = if region_length > 0.0 && cohort.num_samples > 0 {
                (total_bases as f64) / (region_length * (cohort.num_samples as f64))
            } else {
                0.0
            };

            // Clear memory cache to keep memory usage low during batch processing
            read_provider.clear_cache();

            // Store path for dictionary (relative by default, or with prefix)
            // Fix: We want to point to the .lnc_cache.bin file in dictionary.json
            let filename = cache_path.file_name()
                .map(|f| f.to_string_lossy().to_string())
                .unwrap_or_else(|| {
                    let mut p = PathBuf::from(output_bam_path);
                    p.set_extension("lnc_cache.bin");
                    p.file_name().map(|f| f.to_string_lossy().to_string()).unwrap_or_else(|| output_bam_path.clone())
                });
            
            let stored_path = if let Some(ref prefix) = config.path_prefix {
                format!("{}{}", prefix, filename)
            } else {
                filename
            };

            let mut cohort_info = serde_json::Map::new();
            cohort_info.insert("path".to_string(), json!(stored_path));
            if let Some(ref t) = cohort.tissue {
                cohort_info.insert("tissue".to_string(), json!(t));
            }
            if let Some(ref s) = cohort.status {
                cohort_info.insert("status".to_string(), json!(s));
            }
            cohort_info.insert("num_samples".to_string(), json!(cohort.num_samples));
            cohort_info.insert("num_reads".to_string(), json!(num_reads));
            cohort_info.insert("cache_size_bytes".to_string(), json!(cache_size));
            cohort_info.insert("avg_coverage_per_sample".to_string(), json!(avg_coverage_per_sample));

            region_to_cohort_map.entry(key.clone()).or_insert_with(serde_json::Map::new)
                .insert(cohort.name.clone(), json!(cohort_info));
        }

        println!("\n--- Cohort Summary: {} ---", cohort.name);
        println!("  |-- Total Regions:   {}", dict_metadata.len());
        println!("  |-- Elapsed Time:    {:?}", t_cohort_begin.elapsed());
        
        read_provider.benchmark.record_stage(&format!("Cohort: {}", cohort.name), t_cohort_begin.elapsed());
        read_provider.benchmark.set_current_cohort(None);
    }

    // Generate global dictionary.json if output_dir is specified
    if let Some(ref dir) = config.output_dir {
        let mut dict = serde_json::Map::new();
        
        for (key, offset, region_str) in dict_metadata {
            let mut entry = serde_json::Map::new();
            let mut metadata = serde_json::Map::new();
            metadata.insert("offset".to_string(), json!(offset));
            metadata.insert("region".to_string(), json!(region_str));
            if config.compress_headers {
                metadata.insert("experimental_header_compression".to_string(), json!(true));
            }
            entry.insert("metadata".to_string(), json!(metadata));
            
            // Get the cohorts specifically for this region
            let cohort_map = region_to_cohort_map.get(&key).cloned().unwrap_or_else(serde_json::Map::new);
            entry.insert("cohorts".to_string(), json!(cohort_map));
            
            dict.insert(key, json!(entry));
        }

        let dict_path = Path::new(dir).join("dictionary.json");
        let dict_json = serde_json::to_string_pretty(&dict)?;
        fs::write(&dict_path, dict_json).map_err(|e| format!("Failed to write dictionary.json: {}", e))?;
        println!("\nGlobal dictionary written to {}", dict_path.display());
    }

    // Print final benchmarking report if anything was recorded
    let benchmark = lnc_seeker_lib::reads_manager::get_read_provider().benchmark.clone();
    benchmark.print_report();
    let _ = benchmark.write_csv("benchmarking/data/benchmark_report.csv");

    Ok(())
}
