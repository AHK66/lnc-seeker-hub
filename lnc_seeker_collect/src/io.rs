// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use lnc_seeker_lib::pipeline_config::AnyResult;
use noodles::sam;
use noodles::sam::alignment::io::Write;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::Map;

/// Sort records by `(reference_sequence_id, alignment_start)` so they are
/// suitable for downstream BAM indexing (BAI/CSI builders expect
/// coordinate/order grouped by reference id).
pub fn sort_records_for_indexing(records: &mut [RecordBuf]) {
    records.sort_by_key(|r| {
        let rid = r.reference_sequence_id().unwrap_or(usize::MAX);
        let pos = r
            .alignment_start()
            .map(|p| usize::from(p))
            .unwrap_or(usize::MAX);
        (rid, pos)
    });
}

/// Write `records` to `out_bam` using `out_header`, performing basic
/// validation of reference IDs, clearing invalid mate refs, and creating
/// a BAI index. Returns the number of records successfully written.
pub fn write_bam_and_bai(
    out_bam: &str,
    out_header: &sam::Header,
    records: &[RecordBuf],
) -> AnyResult<usize> {
    let out_ref_count = out_header.reference_sequences().len();

    // Clone and ensure required @RG entries exist for intron-based marking.
    let mut header = out_header.clone();

    // Add three read-groups if they're not already present. These IDs are
    // referenced by the `RG` auxiliary tag inserted into records by the
    // intron extraction helpers.
    let rg_ids = ["INTRON_SAME", "INTRON_DIFF", "INTRON_PARTIAL"];
    for id in rg_ids.iter() {
        if !header.read_groups().contains_key(&id.as_bytes()[..]) {
            header.read_groups_mut().insert(String::from(*id).into(), Map::<ReadGroup>::default());
        }
    }

    let out_file = std::fs::File::create(out_bam)
        .map_err(|e| format!("{}: {}", out_bam, e))?;
    let mut writer = noodles::bam::io::Writer::new(noodles::bgzf::io::Writer::new(out_file));
    writer
        .write_header(&header)
        .map_err(|e| format!("{}: {}", out_bam, e))?;

    let mut written = 0usize;
    for record in records.iter() {
        if let Some(rid) = record.reference_sequence_id() {
            if rid >= out_ref_count {
                let q = record
                    .name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                    .unwrap_or_else(|| "<unnamed>".to_string());
                eprintln!(
                    "Skipping write of {}: reference id {} out of range for output header",
                    q, rid
                );
                continue;
            }
        }

        if let Some(mrid) = record.mate_reference_sequence_id() {
            if mrid >= out_ref_count {
                let q = record
                    .name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                    .unwrap_or_else(|| "<unnamed>".to_string());
                eprintln!(
                    "Clearing mate ref for {}: mate ref id {} out of range for output header",
                    q, mrid
                );

                let mut cloned = record.clone();
                *cloned.mate_reference_sequence_id_mut() = None;
                writer
                    .write_alignment_record(&header, &cloned)
                    .map_err(|e| format!("{}: {}", out_bam, e))?;
                written += 1;
                continue;
            }
        }

        writer
            .write_alignment_record(&header, record)
            .map_err(|e| format!("{}: {}", out_bam, e))?;
        written += 1;
    }

    writer
        .try_finish()
        .map_err(|e| format!("{}: {}", out_bam, e))?;

    let index = noodles::bam::fs::index(out_bam)
        .map_err(|e| format!("{}: {}", out_bam, e))?;
    let bai_path = format!("{}.bai", out_bam);
    noodles::bam::bai::fs::write(&bai_path, &index)
        .map_err(|e| format!("{}: {}", bai_path, e))?;

    println!("Wrote {} and index {} ({} reads)", out_bam, bai_path, written);

    Ok(written)
}
