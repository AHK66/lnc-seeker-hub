// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::collections::HashMap;
use std::io;
use serde::{Deserialize, Serialize};
use noodles::sam;
use sam::alignment::Record;
use sam::alignment::record::data::field::Tag;
use sam::alignment::record::data::field::Value;

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct JunctionSpan {
    pub reference: String,
    pub start: usize,
    pub end: usize,
    pub reads: usize,
    pub reads_clean: usize,
    pub anchored_start: usize,
    pub anchored_end: usize,
    #[serde(default = "default_junction_strand")]
    pub junction_strand: String,
    #[serde(default = "default_strand_source")]
    pub strand_source: String,
}

fn default_junction_strand() -> String {
    ".".to_string()
}

fn default_strand_source() -> String {
    "fallback".to_string()
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct ReadSegment {
    pub start: usize,
    pub end: usize,
    pub is_mate: bool,
    pub mismatches: u8,
    pub insertions: u8,
    pub is_followed_by_deletion: bool,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct ReadInfo {
    pub name: String,
    pub reference: String,
    pub mapping_quality: u8,
    #[serde(default)]
    pub strand: i8,
    pub start: usize,
    pub end: usize,
    pub segments: Vec<ReadSegment>,
}

impl ReadInfo {
    fn strand_char_from_value(value: Value<'_>) -> Option<char> {
        match value {
            Value::Character(b) => {
                let c = char::from(b);
                if c == '+' || c == '-' { Some(c) } else { None }
            }
            Value::String(s) => {
                let c = s.iter().next().copied().map(char::from).unwrap_or('.');
                if c == '+' || c == '-' { Some(c) } else { None }
            }
            _ => None,
        }
    }

    fn infer_splice_strand<R>(record: &R) -> i8
    where
        R: Record,
    {
        let data = record.data();

        if let Some(Ok(v)) = data.get(&Tag::new(b'X', b'S')) {
            if let Some(c) = Self::strand_char_from_value(v) {
                return if c == '+' { 1 } else { -1 };
            }
        }

        if let Some(Ok(v)) = data.get(&Tag::new(b't', b's')) {
            if let Some(c) = Self::strand_char_from_value(v) {
                if let Ok(flags) = record.flags() {
                    let reverse = flags.is_reverse_complemented();
                    if reverse {
                        return if c == '+' { -1 } else { 1 };
                    }
                }
                return if c == '+' { 1 } else { -1 };
            }
        }

        0
    }

    pub fn from_record<R>(record: &R, header: &sam::Header, ref_seq: Option<&[u8]>) -> io::Result<Option<Self>>
    where
        R: Record,
    {
        let Some(ref_id) = record.reference_sequence_id(header).transpose()? else {
            return Ok(None);
        };

        let (ref_name, _) = header.reference_sequences().get_index(ref_id).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "Invalid reference sequence ID")
        })?;

        let Some(alignment_start) = record.alignment_start().transpose()? else {
            return Ok(None);
        };

        let mapping_quality = record.mapping_quality().transpose()?.map(u8::from).unwrap_or(60);
        let read_name = record.name().map(|n| n.to_string()).unwrap_or_else(|| "unknown".to_string());
        let is_mate = record.flags()?.is_last_segment();
        let strand = Self::infer_splice_strand(record);
        
        let mut curr = usize::from(alignment_start) - 1;
        let mut read_pos = 0;
        let mut segments = Vec::new();
        let mut current_segment_start = Some(curr);
        let mut current_mismatches = 0u8;
        let mut current_insertions = 0u8;

        let seq = record.sequence();

        for op_result in record.cigar().iter() {
            let op = op_result?;
            let len = op.len();
            match op.kind() {
                sam::alignment::record::cigar::op::Kind::Match
                | sam::alignment::record::cigar::op::Kind::SequenceMatch => {
                    if current_segment_start.is_none() {
                        current_segment_start = Some(curr);
                    }
                    if let Some(rs) = ref_seq {
                        for i in 0..len {
                            let r_pos = read_pos + i;
                            let g_pos = curr + i;
                            if r_pos < seq.len() && g_pos < rs.len() {
                                let read_base = seq.get(r_pos).map(u8::from).unwrap_or(0);
                                let ref_base = rs[g_pos];
                                // Basic comparison (case insensitive for genome)
                                if (read_base != 0) && (read_base.to_ascii_uppercase() != ref_base.to_ascii_uppercase()) {
                                    current_mismatches = current_mismatches.saturating_add(1);
                                }
                            }
                        }
                    }
                    curr += len;
                    read_pos += len;
                }
                sam::alignment::record::cigar::op::Kind::SequenceMismatch => {
                    if current_segment_start.is_none() {
                        current_segment_start = Some(curr);
                    }
                    current_mismatches = current_mismatches.saturating_add(len as u8);
                    curr += len;
                    read_pos += len;
                }
                sam::alignment::record::cigar::op::Kind::Insertion => {
                    current_insertions = current_insertions.saturating_add(len as u8).min(7);
                    read_pos += len;
                }
                sam::alignment::record::cigar::op::Kind::Deletion => {
                    if let Some(s) = current_segment_start {
                        segments.push(ReadSegment { 
                            start: s, 
                            end: curr, 
                            is_mate,
                            mismatches: current_mismatches.min(15),
                            insertions: current_insertions,
                            is_followed_by_deletion: true,
                        });
                        current_segment_start = None;
                        current_mismatches = 0;
                        current_insertions = 0;
                    }
                    curr += len;
                }
                sam::alignment::record::cigar::op::Kind::Skip => {
                    if let Some(s) = current_segment_start {
                        segments.push(ReadSegment { 
                            start: s, 
                            end: curr, 
                            is_mate,
                            mismatches: current_mismatches.min(15),
                            insertions: current_insertions,
                            is_followed_by_deletion: false,
                        });
                        current_segment_start = None;
                        current_mismatches = 0;
                        current_insertions = 0;
                    }
                    curr += len;
                }
                sam::alignment::record::cigar::op::Kind::SoftClip => {
                    read_pos += len;
                }
                sam::alignment::record::cigar::op::Kind::HardClip => {}
                _ => {}
            }
        }
        
        if let Some(s) = current_segment_start {
            segments.push(ReadSegment { 
                start: s, 
                end: curr, 
                is_mate,
                mismatches: current_mismatches,
                insertions: current_insertions,
                is_followed_by_deletion: false,
            });
        }

        let alignment_end = record.alignment_end().transpose()?.map(usize::from).unwrap_or(curr);
        let alignment_start_0 = usize::from(alignment_start) - 1;

        Ok(Some(ReadInfo {
            name: read_name,
            reference: ref_name.to_string(),
            mapping_quality,
            strand,
            start: alignment_start_0,
            end: alignment_end,
            segments,
        }))
    }
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct JunctionPoint {
    pub reference: String,
    pub position: usize,
    pub junctions: usize,
    pub reads: usize,
    pub avg_before: f64,
    pub avg_after: f64,
    pub change_pct: f64,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Annotation {
    pub reference: String,
    pub start: usize,
    pub end: usize,
    pub feature: String,
    pub gene_id: String,
    pub gene_name: Option<String>,
    pub transcript_id: String,
    pub exon_number: Option<String>,
    pub strand: String,
}

/// Binary-optimized version of Annotation for memory-mapped storage and O(log n) lookups.
/// Uses indices into a global string pool to eliminate redundant allocations.
#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
pub struct PackedAnnotation {
    pub chrom_idx: u16,
    pub start: u32,
    pub end: u32,
    pub feature_idx: u8,       // 0=unknown, 1=exon, 2=transcript, 3=cds, 4=gene, 5=UTR
    pub strand: u8,            // 0=+, 1=-, 2=.
    pub gene_id_idx: u32,
    pub gene_name_idx: u32,    // u32::MAX for None
    pub transcript_id_idx: u32,
    pub exon_number_idx: u32,  // u32::MAX for None
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct LbaFile {
    pub magic: [u8; 4],        // b"LBA\x01"
    pub chroms: Vec<String>,
    pub string_pool: Vec<String>,
    pub records: Vec<PackedAnnotation>,
}

impl LbaFile {
    pub fn estimate_size_bytes(&self) -> usize {
        let records_size = self.records.len() * std::mem::size_of::<PackedAnnotation>();
        let chroms_size: usize = self.chroms.iter().map(|s| s.len() + 24).sum();
        let pool_size: usize = self.string_pool.iter().map(|s| s.len() + 24).sum();
        records_size + chroms_size + pool_size + 64 // 64 bytes for fixed overhead
    }
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct SampleResult {
    pub reference: String,
    pub positions: Vec<usize>,
    pub depths: Vec<u32>,
    pub depths_hq: Vec<u32>,
    pub depths_ambiguity: Vec<u32>,
    pub junction_points: Vec<JunctionPoint>,
    pub junction_spans: Vec<JunctionSpan>,
    #[serde(default)]
    pub min_x: usize,
    #[serde(default)]
    pub max_x: usize,
    #[serde(default)]
    pub coverage_memory_mb: f64,
}

impl SampleResult {
    /// Detects the core area of the plot (excluding outliers) and sets min_x/max_x.
    /// Does not trim the data vectors.
    pub fn detect_outliers(&mut self) {
        if self.depths.is_empty() {
            return;
        }

        // Use depths (gray line) to determine the main bulk of data.
        let total_depth: u64 = self.depths.iter().map(|&d| d as u64).sum();
        if total_depth == 0 {
            self.min_x = self.positions[0];
            self.max_x = *self.positions.last().unwrap_or(&self.positions[0]);
            return;
        }

        // We want to identify the core 99.99% of the data.
        let threshold = (total_depth as f64 * 0.0001) as u64;

        let mut current_sum = 0u64;
        let mut min_idx = 0;
        for (i, &depth) in self.depths.iter().enumerate() {
            current_sum += depth as u64;
            if current_sum >= threshold {
                min_idx = i;
                break;
            }
        }

        current_sum = 0;
        let mut max_idx = self.depths.len() - 1;
        for (i, &depth) in self.depths.iter().enumerate().rev() {
            current_sum += depth as u64;
            if current_sum >= threshold {
                max_idx = i;
                break;
            }
        }

        if min_idx >= max_idx {
            self.min_x = self.positions[0];
            self.max_x = *self.positions.last().unwrap_or(&self.positions[0]);
            return;
        }

        self.min_x = self.positions[min_idx];
        self.max_x = self.positions[max_idx];
    }

    pub fn filter_outliers(&mut self) {
        self.detect_outliers();
        
        let start_pos = self.min_x;
        let end_pos = self.max_x;

        // Find indices for trimming
        let min_idx = self.positions.iter().position(|&p| p == start_pos).unwrap_or(0);
        let max_idx = self.positions.iter().rposition(|&p| p == end_pos).unwrap_or(self.positions.len() - 1);

        // Trim primary vectors
        self.positions = self.positions[min_idx..=max_idx].to_vec();
        self.depths = self.depths[min_idx..=max_idx].to_vec();
        self.depths_hq = self.depths_hq[min_idx..=max_idx].to_vec();
        self.depths_ambiguity = self.depths_ambiguity[min_idx..=max_idx].to_vec();

        // Filter junctions
        self.junction_points.retain(|jp| jp.position >= start_pos && jp.position <= end_pos);
        self.junction_spans.retain(|js| js.start >= start_pos && js.end <= end_pos);
    }
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct AnalysisResult {
    pub samples: HashMap<String, SampleResult>,
    pub annotations: Vec<Annotation>,
    #[serde(default)]
    pub min_x: usize,
    #[serde(default)]
    pub max_x: usize,
    #[serde(default)]
    pub cache_core_mb: f64,
    #[serde(default)]
    pub cache_related_mb: f64,
    #[serde(default)]
    pub cache_annotation_mb: f64,
}
