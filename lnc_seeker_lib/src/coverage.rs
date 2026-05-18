// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::collections::{BTreeMap, HashMap};
use std::io;
use noodles::{bam, sam};
use sam::alignment::Record as _;
use sam::alignment::record::data::field::Tag;
use sam::alignment::record::data::field::Value;
use crate::models::{JunctionPoint, JunctionSpan};
use crate::reads_manager::{CompactRead, CompactSegment};

fn infer_splice_strand_from_bam_record(record: &bam::Record) -> i8 {
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

    let data = record.data();

    if let Some(Ok(v)) = data.get(&Tag::new(b'X', b'S')) {
        if let Some(c) = strand_char_from_value(v) {
            return if c == '+' { 1 } else { -1 };
        }
    }

    if let Some(Ok(v)) = data.get(&Tag::new(b't', b's')) {
        if let Some(c) = strand_char_from_value(v) {
            let reverse = record.flags().is_reverse_complemented();
            if reverse {
                return if c == '+' { -1 } else { 1 };
            }
            return if c == '+' { 1 } else { -1 };
        }
    }

    0
}

pub struct JunctionStore {
    pub data: HashMap<(usize, usize, usize), (usize, usize, usize, usize, usize)>,
    pub start_anchors: HashMap<(usize, usize), usize>,
    pub end_anchors: HashMap<(usize, usize), usize>,
}

impl Default for JunctionStore {
    fn default() -> Self {
        Self {
            data: HashMap::new(),
            start_anchors: HashMap::new(),
            end_anchors: HashMap::new(),
        }
    }
}

impl JunctionStore {
    pub fn add_junction(&mut self, ref_id: usize, start: usize, end: usize, is_clean: bool, strand: i8) {
        let entry = self.data.entry((ref_id, start, end)).or_insert((0, 0, 0, 0, 0));
        entry.0 += 1;
        if is_clean {
            entry.1 += 1;
        }
        match strand {
            1 => entry.2 += 1,
            -1 => entry.3 += 1,
            _ => entry.4 += 1,
        }
    }

    pub fn add_anchor(&mut self, ref_id: usize, start: usize, end: usize) {
        *self.start_anchors.entry((ref_id, end)).or_default() += 1;
        *self.end_anchors.entry((ref_id, start)).or_default() += 1;
    }

    pub fn get_points(&self, ref_names: &[String], coverage: &Coverage, hq: bool) -> Vec<JunctionPoint> {
        let mut results = Vec::new();
        let mut ref_points: BTreeMap<usize, BTreeMap<usize, (usize, usize)>> = BTreeMap::new();
        for ((ref_id, start, end), (count, _, _, _, _)) in &self.data {
            let points = ref_points.entry(*ref_id).or_default();
            let s_entry = points.entry(*start).or_insert((0, 0));
            s_entry.0 += 1;
            s_entry.1 += *count;
            let e_entry = points.entry(*end).or_insert((0, 0));
            e_entry.0 += 1;
            e_entry.1 += *count;
        }

        for (ref_id, points) in ref_points {
            let name = ref_names
                .get(ref_id)
                .cloned()
                .unwrap_or_else(|| format!("{}", ref_id));

            for (pos, (j_count, r_count)) in points {
                let avg_before = coverage.get_average_coverage(ref_id, pos.saturating_sub(5), pos, hq);
                let avg_after = coverage.get_average_coverage(ref_id, pos, pos + 5, hq);
                let change = if avg_before > 0.0 {
                    (avg_after - avg_before) / avg_before
                } else if avg_after > 0.0 {
                    10.0 // Cap at 1000% if before was 0
                } else {
                    0.0
                };

                results.push(JunctionPoint {
                    reference: name.clone(),
                    position: pos,
                    junctions: j_count,
                    reads: r_count,
                    avg_before,
                    avg_after,
                    change_pct: change * 100.0,
                });
            }
        }
        results
    }

    pub fn get_spans(&self, ref_names: &[String]) -> Vec<JunctionSpan> {
        let mut results = Vec::new();
        for ((ref_id, start, end), (count, count_clean, count_plus, count_minus, _count_unknown)) in &self.data {
            let name = ref_names
                .get(*ref_id)
                .cloned()
                .unwrap_or_else(|| format!("{}", ref_id));

            let anchored_start = self.start_anchors.get(&(*ref_id, *start)).copied().unwrap_or(0);
            let anchored_end = self.end_anchors.get(&(*ref_id, *end)).copied().unwrap_or(0);

            let junction_strand = if *count_plus > *count_minus {
                "+"
            } else if *count_minus > *count_plus {
                "-"
            } else {
                "."
            };

            let strand_source = if junction_strand == "." {
                "fallback"
            } else {
                "tag"
            };

            results.push(JunctionSpan {
                reference: name,
                start: *start,
                end: *end,
                reads: *count,
                reads_clean: *count_clean,
                anchored_start,
                anchored_end,
                junction_strand: junction_strand.to_string(),
                strand_source: strand_source.to_string(),
            });
        }
        results
    }
}

#[derive(Clone)]
pub struct ReferenceCoverage {
    pub offset: usize,
    pub counts: Vec<u32>,
}

impl ReferenceCoverage {
    pub fn size_bytes(&self) -> usize {
        std::mem::size_of::<Self>() + (self.counts.capacity() * 4)
    }

    pub fn new(start: usize, end: usize) -> Self {
        Self {
            offset: start,
            counts: vec![0; end - start],
        }
    }

    pub fn expand(&mut self, start: usize, end: usize) {
        if start < self.offset {
            let padding = self.offset - start;
            let mut new_counts = vec![0; padding];
            new_counts.extend(self.counts.drain(..));
            self.counts = new_counts;
            self.offset = start;
        }

        let current_end = self.offset + self.counts.len();
        if end > current_end {
            let padding = end - current_end;
            self.counts.extend(std::iter::repeat(0).take(padding));
        }
    }

    pub fn increment_range(&mut self, start: usize, end: usize) {
        let relative_start = start - self.offset;
        let relative_end = end - self.offset;
        for i in relative_start..relative_end {
            self.counts[i] += 1;
        }
    }

    pub fn get(&self, position_0: usize) -> u32 {
        if position_0 >= self.offset {
            let relative_pos = position_0 - self.offset;
            self.counts.get(relative_pos).copied().unwrap_or(0)
        } else {
            0
        }
    }
}

pub struct Coverage {
    pub data: Vec<Option<ReferenceCoverage>>,
    pub data_hq: Vec<Option<ReferenceCoverage>>,
    pub data_ambiguity: Vec<Option<ReferenceCoverage>>,
}

impl Coverage {
    pub fn size_bytes(&self) -> usize {
        let mut total = 0;
        for c in self.data.iter().flatten() {
            total += c.size_bytes();
        }
        for c in self.data_hq.iter().flatten() {
            total += c.size_bytes();
        }
        for c in self.data_ambiguity.iter().flatten() {
            total += c.size_bytes();
        }
        total
    }

    pub fn new(num_references: usize) -> Self {
        let data = vec![None; num_references];
        let data_hq = vec![None; num_references];
        let data_ambiguity = vec![None; num_references];
        Self { data, data_hq, data_ambiguity }
    }

    pub fn add_record(&mut self, record: &bam::Record, junctions: &mut JunctionStore, junctions_hq: &mut JunctionStore, min_mapping_quality: u8, ambiguity_min_mapping_quality: u8) -> io::Result<()> {
        let Some(reference_sequence_id) = record.reference_sequence_id().transpose()? else {
            return Ok(());
        };

        // Standard Filter: Skip non-primary, supplementary, and duplicate reads for all coverage tracks.
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() || flags.is_duplicate() {
            return Ok(());
        }

        let Some(alignment_start) = record.alignment_start().transpose()? else {
            return Ok(());
        };

        let Some(alignment_end) = record.alignment_end().transpose()? else {
            return Ok(());
        };

        let start = usize::from(alignment_start) - 1;
        let end = usize::from(alignment_end);

        let mapping_quality = record.mapping_quality().map(u8::from).unwrap_or(60); // Default high if missing
        let is_high_quality = mapping_quality >= min_mapping_quality;
        let is_ambiguity_quality = mapping_quality >= ambiguity_min_mapping_quality;
        let read_strand = infer_splice_strand_from_bam_record(record);

        if let Some(reference_slot) = self.data.get_mut(reference_sequence_id) {
            if let Some(rc) = reference_slot {
                rc.expand(start, end);
            } else {
                *reference_slot = Some(ReferenceCoverage::new(start, end));
            }
        }

        if is_high_quality {
            if let Some(reference_slot_hq) = self.data_hq.get_mut(reference_sequence_id) {
                if let Some(rc_hq) = reference_slot_hq {
                    rc_hq.expand(start, end);
                } else {
                    *reference_slot_hq = Some(ReferenceCoverage::new(start, end));
                }
            }
        }

        if is_ambiguity_quality {
            if let Some(reference_slot_amb) = self.data_ambiguity.get_mut(reference_sequence_id) {
                if let Some(rc_amb) = reference_slot_amb {
                    rc_amb.expand(start, end);
                } else {
                    *reference_slot_amb = Some(ReferenceCoverage::new(start, end));
                }
            }
        }

        junctions.add_anchor(reference_sequence_id, start, end);
        if is_high_quality {
            junctions_hq.add_anchor(reference_sequence_id, start, end);
        }

        let mut curr = start;
        for result in record.cigar().iter() {
            let op = result?;
            let len = op.len();
            match op.kind() {
                sam::alignment::record::cigar::op::Kind::Match
                | sam::alignment::record::cigar::op::Kind::SequenceMatch
                | sam::alignment::record::cigar::op::Kind::SequenceMismatch => {
                    if let Some(Some(rc)) = self.data.get_mut(reference_sequence_id) {
                        rc.increment_range(curr, curr + len);
                    }
                    if is_high_quality {
                        if let Some(Some(rc_hq)) = self.data_hq.get_mut(reference_sequence_id) {
                            rc_hq.increment_range(curr, curr + len);
                        }
                    }
                    if is_ambiguity_quality {
                        if let Some(Some(rc_amb)) = self.data_ambiguity.get_mut(reference_sequence_id) {
                            rc_amb.increment_range(curr, curr + len);
                        }
                    }
                    curr += len;
                }
                sam::alignment::record::cigar::op::Kind::Deletion => {
                    curr += len;
                }
                sam::alignment::record::cigar::op::Kind::Skip => {
                    junctions.add_junction(reference_sequence_id, curr, curr + len, true, read_strand);
                    if is_high_quality {
                        junctions_hq.add_junction(reference_sequence_id, curr, curr + len, true, read_strand);
                    }
                    curr += len;
                }
                _ => {}
            }
        }

        Ok(())
    }

    pub fn add_compact_read(
        &mut self,
        read: &CompactRead,
        segments: &[CompactSegment],
        segment_tags: &[u8],
        junctions: &mut JunctionStore,
        junctions_hq: &mut JunctionStore,
        min_mapping_quality: u8,
        ambiguity_min_mapping_quality: u8,
    ) {
        let reference_sequence_id = read.ref_id as usize;
        let start = read.start as usize;
        let end = read.end as usize;

        let mapping_quality = read.mq;
        let is_high_quality = mapping_quality >= min_mapping_quality;
        let is_ambiguity_quality = mapping_quality >= ambiguity_min_mapping_quality;
        let read_strand = match (read.flags >> 2) & 0x03 {
            1 => 1,
            2 => -1,
            _ => 0,
        };

        if let Some(reference_slot) = self.data.get_mut(reference_sequence_id) {
            if let Some(rc) = reference_slot {
                rc.expand(start, end);
            } else {
                *reference_slot = Some(ReferenceCoverage::new(start, end));
            }
        }

        if is_high_quality {
            if let Some(reference_slot_hq) = self.data_hq.get_mut(reference_sequence_id) {
                if let Some(rc_hq) = reference_slot_hq {
                    rc_hq.expand(start, end);
                } else {
                    *reference_slot_hq = Some(ReferenceCoverage::new(start, end));
                }
            }
        }

        if is_ambiguity_quality {
            if let Some(reference_slot_amb) = self.data_ambiguity.get_mut(reference_sequence_id) {
                if let Some(rc_amb) = reference_slot_amb {
                    rc_amb.expand(start, end);
                } else {
                    *reference_slot_amb = Some(ReferenceCoverage::new(start, end));
                }
            }
        }

        junctions.add_anchor(reference_sequence_id, start, end);
        if is_high_quality {
            junctions_hq.add_anchor(reference_sequence_id, start, end);
        }

        for i in 0..read.seg_count as usize {
            let seg = &segments[read.seg_offset as usize + i];
            let seg_start = seg.start as usize;
            let seg_end = seg.end as usize;

            if let Some(Some(rc)) = self.data.get_mut(reference_sequence_id) {
                rc.increment_range(seg_start, seg_end);
            }
            if is_high_quality {
                if let Some(Some(rc_hq)) = self.data_hq.get_mut(reference_sequence_id) {
                    rc_hq.increment_range(seg_start, seg_end);
                }
            }
            if is_ambiguity_quality {
                if let Some(Some(rc_amb)) = self.data_ambiguity.get_mut(reference_sequence_id) {
                    rc_amb.increment_range(seg_start, seg_end);
                }
            }

            if i < (read.seg_count as usize).saturating_sub(1) {
                let next_seg = &segments[read.seg_offset as usize + i + 1];
                let tag = segment_tags[read.seg_offset as usize + i];
                let next_tag = segment_tags[read.seg_offset as usize + i + 1];
                let is_deletion = (tag & 0x01) != 0;

                if !is_deletion {
                    let j_start = seg_end;
                    let j_end = next_seg.start as usize;
                    let is_clean = (tag >> 1 == 0) && (next_tag >> 1 == 0);
                    
                    junctions.add_junction(reference_sequence_id, j_start, j_end, is_clean, read_strand);
                    if is_high_quality {
                        junctions_hq.add_junction(reference_sequence_id, j_start, j_end, is_clean, read_strand);
                    }
                }
            }
        }
    }

    pub fn get_average_coverage(&self, reference_sequence_id: usize, start: usize, end: usize, hq: bool) -> f64 {
        if start >= end {
            return 0.0;
        }
        let data = if hq { &self.data_hq } else { &self.data };
        let Some(Some(rc)) = data.get(reference_sequence_id) else {
            return 0.0;
        };

        let mut sum = 0u64;
        let count = end - start;
        for pos in start..end {
            sum += rc.get(pos) as u64;
        }
        sum as f64 / count as f64
    }

    pub fn get_points(&self, ref_names: &[String]) -> (String, Vec<usize>, Vec<u32>, Vec<u32>, Vec<u32>) {
        let mut best_ref_id = None;
        let mut max_depth = 0u64;

        for (i, rc_opt) in self.data.iter().enumerate() {
            if let Some(rc) = rc_opt {
                let d: u64 = rc.counts.iter().map(|&x| x as u64).sum();
                if d > max_depth {
                    max_depth = d;
                    best_ref_id = Some(i);
                }
            }
        }

        let mut positions = Vec::new();
        let mut depths = Vec::new();
        let mut depths_hq = Vec::new();
        let mut depths_ambiguity = Vec::new();
        let mut ref_name = String::from("unknown");

        if let Some(i) = best_ref_id {
            ref_name = ref_names
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("{}", i));

            if let Some(rc) = &self.data[i] {
                let rc_hq_opt = self.data_hq.get(i).and_then(|opt| opt.as_ref());
                let rc_amb_opt = self.data_ambiguity.get(i).and_then(|opt| opt.as_ref());

                for (offset, &depth) in rc.counts.iter().enumerate() {
                    let pos = rc.offset + offset;
                    let depth_hq = rc_hq_opt.map(|rc_hq| rc_hq.get(pos)).unwrap_or(0);
                    let depth_amb = rc_amb_opt.map(|rc_amb| rc_amb.get(pos)).unwrap_or(0);
                    positions.push(pos);
                    depths.push(depth);
                    depths_hq.push(depth_hq);
                    depths_ambiguity.push(depth_amb);
                }
            }
        }

        (ref_name, positions, depths, depths_hq, depths_ambiguity)
    }

    pub fn finalize(
        &self, 
        ref_names: &[String],
        _junctions: &JunctionStore,
        junctions_hq: &JunctionStore,
        sample_res: &mut crate::models::SampleResult
    ) -> std::collections::HashMap<String, (usize, usize)> {
        let (main_ref, positions, depths, depths_hq, depths_ambiguity) = self.get_points(ref_names);
        
        sample_res.reference = main_ref.clone();
        sample_res.positions = positions;
        sample_res.depths = depths;
        sample_res.depths_hq = depths_hq;
        sample_res.depths_ambiguity = depths_ambiguity;
        sample_res.junction_points = junctions_hq.get_points(ref_names, self, true);
        sample_res.junction_spans = junctions_hq.get_spans(ref_names);
        
        // Filter junctions to main reference
        sample_res.junction_points.retain(|jp| jp.reference == main_ref);
        sample_res.junction_spans.retain(|js| js.reference == main_ref);

        let mut covered_ranges = std::collections::HashMap::new();
        for (i, rc_opt) in self.data.iter().enumerate() {
            if let Some(rc) = rc_opt {
                let name = ref_names.get(i)
                    .cloned()
                    .unwrap_or_else(|| format!("{}", i));
                covered_ranges.insert(name, (rc.offset, rc.offset + rc.counts.len()));
            }
        }
        covered_ranges
    }
}
