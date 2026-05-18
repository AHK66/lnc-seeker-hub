// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
#[cfg(test)]
mod tests {
    use crate::analysis::run_analysis;
    use crate::config::Config;
    use crate::models::{ReadInfo, ReadSegment};
    use crate::processing::get_junction_reads;
    use crate::reads_manager::{get_read_provider, CacheStatus};
    use noodles::bam;
    use std::fs::File;
    use std::io::BufReader;
    use std::sync::Arc;

    #[test]
    fn test_cache_integrity() {
        let provider = get_read_provider();
        provider.clear_cache();

        let mut reads = Vec::new();
        let ref_name = "chr1".to_string();

        reads.push(ReadInfo {
            name: "read_A".to_string(),
            reference: ref_name.clone(),
            mapping_quality: 60,
            strand: 0,
            start: 50,
            end: 250,
            segments: vec![
                ReadSegment {
                    start: 50,
                    end: 100,
                    is_mate: false,
                    mismatches: 0,
                    insertions: 0,
                    is_followed_by_deletion: false,
                },
                ReadSegment {
                    start: 200,
                    end: 250,
                    is_mate: false,
                    mismatches: 0,
                    insertions: 0,
                    is_followed_by_deletion: false,
                },
            ],
        });

        reads.push(ReadInfo {
            name: "read_B".to_string(),
            reference: ref_name.clone(),
            mapping_quality: 60,
            strand: 0,
            start: 80,
            end: 100,
            segments: vec![ReadSegment {
                start: 80,
                end: 100,
                is_mate: false,
                mismatches: 0,
                insertions: 0,
                is_followed_by_deletion: false,
            }],
        });

        reads.push(ReadInfo {
            name: "read_B".to_string(),
            reference: ref_name.clone(),
            mapping_quality: 60,
            strand: 0,
            start: 200,
            end: 220,
            segments: vec![ReadSegment {
                start: 200,
                end: 220,
                is_mate: true,
                mismatches: 0,
                insertions: 0,
                is_followed_by_deletion: false,
            }],
        });

        let bam_path = "test_mock.bam";
        let ref_names = vec![ref_name.clone()];
        let ref_lengths = vec![1000usize];
        let _ = provider.commit_caching_compact(
            bam_path,
            reads,
            ref_names,
            ref_lengths,
            CacheStatus::Complete,
        );

        let results = provider.get_filtered_reads(bam_path, &ref_name, 100, 200, 20);
        assert!(results.is_some());
        let (filtered, status) = results.unwrap();
        assert_eq!(status, CacheStatus::Complete);
        assert_eq!(filtered.len(), 3);

        let names: Vec<String> = filtered.iter().map(|r| r.name.clone()).collect();
        assert!(names.contains(&"read_A".to_string()));
        assert!(names.contains(&"read_B".to_string()));
    }

    #[test]
    #[ignore]
    fn test_real_bam_cache() {
        let bam_path = match std::env::var("LNC_SEEKER_TEST_BAM") {
            Ok(path) => path,
            Err(_) => {
                println!(
                    "Skipping real BAM test: set LNC_SEEKER_TEST_BAM to a BAM path before running"
                );
                return;
            }
        };
        if !std::path::Path::new(&bam_path).exists() {
            println!("Skipping real BAM test: file not found at {}", bam_path);
            return;
        }

        let provider = get_read_provider();
        provider.clear_cache();

        let file = File::open(&bam_path).unwrap();
        let mut reader = bam::io::Reader::new(BufReader::new(file));
        let header = reader.read_header().unwrap();

        let mut target_junction = None;
        for result in reader.records() {
            let record = result.unwrap();
            let info = ReadInfo::from_record(&record, &header, None).unwrap().unwrap();
            if info.segments.len() > 1 {
                target_junction = Some((
                    info.reference.clone(),
                    info.segments[0].end,
                    info.segments[1].start,
                ));
                break;
            }
        }

        let (ref_name, s_target, e_target) = target_junction.expect("No spliced reads found in BAM!");
        println!("Testing junction: {}:{}-{}", ref_name, s_target, e_target);

        let config_json = format!(
            r#"{{
            "data_selection": {{
                "bam_paths": ["{}"],
                "bam_to_cohort": {{}},
                "filter_outliers": true
            }},
            "coverage_and_junctions_profile": {{
                "min_mapping_quality": 20,
                "high_ambiguity_highlighting": {{
                    "ambiguity_min_mapping_quality": 5
                }}
            }}
        }}"#,
            bam_path.replace("\\", "\\\\")
        );

        let config: Config = serde_json::from_str(&config_json).unwrap();
        let progress = Arc::new(crate::progress::ProgressData {
            stage: std::sync::atomic::AtomicU32::new(0),
            current: std::sync::atomic::AtomicU32::new(0),
            total: std::sync::atomic::AtomicU32::new(0),
        });
        run_analysis(&config, progress).unwrap();

        let start_cache = std::time::Instant::now();
        let cache_results =
            get_junction_reads(&bam_path, &ref_name, s_target, e_target, 20, 100, &None, false)
                .unwrap();
        let cache_time = start_cache.elapsed();

        provider.clear_cache();
        let start_disk = std::time::Instant::now();
        let disk_results =
            get_junction_reads(&bam_path, &ref_name, s_target, e_target, 20, 100, &None, false)
                .unwrap();
        let disk_time = start_disk.elapsed();

        println!("Cache hits: {}, Time: {:?}", cache_results.len(), cache_time);
        println!("Disk hits: {}, Time: {:?}", disk_results.len(), disk_time);

        assert_eq!(
            cache_results.len(),
            disk_results.len(),
            "Results length mismatch!"
        );
        assert!(cache_time < disk_time, "Cache should be faster than disk!");
    }
}
