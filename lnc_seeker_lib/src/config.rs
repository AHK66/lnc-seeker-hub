// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::collections::HashMap;
use serde::{Deserialize, Serialize};
use serde::de;
use serde_json::Value as JsonValue;

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct DataSelection {
    #[serde(default)]
    pub bam_path: Option<String>,
    #[serde(default, deserialize_with = "deserialize_bam_paths")]
    pub bam_paths: Vec<String>,
    #[serde(default)]
    pub bam_to_cohort: HashMap<String, String>,
    #[serde(default)]
    pub gtf_path: Option<String>,
    #[serde(default, deserialize_with = "deserialize_bam_paths")]
    pub gtf_paths: Vec<String>,
    #[serde(default)]
    pub genome_path: Option<String>,
    #[serde(default)]
    pub selected_gtfs: Vec<String>,
    pub assembly_report_path: Option<String>,
    #[serde(default = "default_true")]
    pub filter_outliers: bool,
    #[serde(default = "default_true")]
    pub filter_annotations: bool,
    #[serde(default)]
    pub analysis_start: Option<usize>,
    #[serde(default)]
    pub analysis_end: Option<usize>,
    #[serde(default)]
    pub analysis_reference: Option<String>,
    #[serde(default)]
    pub gtf_offsets: HashMap<String, i32>,
    #[serde(default)]
    pub max_cache_memory_mb: Option<f64>,
}

fn deserialize_bam_paths<'de, D>(deserializer: D) -> Result<Vec<String>, D::Error>
where
    D: de::Deserializer<'de>,
{
    let v = JsonValue::deserialize(deserializer).map_err(de::Error::custom)?;
    match v {
        JsonValue::Array(arr) => {
            let mut res = Vec::new();
            for item in arr {
                if let JsonValue::String(s) = item {
                    res.push(s);
                } else {
                    return Err(de::Error::custom("bam_paths array must contain strings"));
                }
            }
            Ok(res)
        }
        JsonValue::Object(map) => {
            let mut res = Vec::new();
            for (_k, val) in map {
                if let JsonValue::Array(arr) = val {
                    for item in arr {
                        if let JsonValue::String(s) = item {
                            res.push(s);
                        } else {
                            return Err(de::Error::custom("bam_paths mapping values must be arrays of strings"));
                        }
                    }
                } else {
                    return Err(de::Error::custom("bam_paths mapping values must be arrays"));
                }
            }
            Ok(res)
        }
        JsonValue::Null => Ok(Vec::new()),
        _ => Err(de::Error::custom("bam_paths must be array of strings or mapping of name -> array")),
    }
}

#[derive(Debug, Deserialize, Serialize, Clone, Default)]
pub struct AmbiguitySettings {
    #[serde(default = "default_min_mapping_quality")]
    pub ambiguity_min_mapping_quality: u8,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct CoverageProfile {
    #[serde(default = "default_min_mapping_quality")]
    pub min_mapping_quality: u8,
    #[serde(rename = "high_ambiguity_highlighting", default)]
    pub ambiguity: AmbiguitySettings,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Config {
    pub data_selection: DataSelection,
    #[serde(rename = "coverage_and_junctions_profile")]
    pub coverage: CoverageProfile,
}

impl Config {
    pub fn get_paths(&self) -> Vec<String> {
        let mut paths = self.data_selection.bam_paths.clone();
        if let Some(path) = &self.data_selection.bam_path {
            if !paths.contains(path) {
                paths.push(path.clone());
            }
        }
        paths
    }
}

pub fn default_min_mapping_quality() -> u8 {
    20
}

pub fn default_true() -> bool {
    true
}
