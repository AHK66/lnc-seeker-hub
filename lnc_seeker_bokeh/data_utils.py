# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import pandas as pd
import math
from lnc_seeker_bokeh.state import log_safe, strip_id

def process_analysis_data(state):
    """Processes the in-memory analysis data from Rust."""
    # Process annotations if we have any data (or an empty structure)
    data = state["analysis_data"] or {}
    
    # Standardize sample names (keys) to cohort names
    samples = data.get("samples", {})
    new_samples = {}
    with state["lock"]:
        stem_map = state.get("stem_to_cohort", {})
        for k, v in samples.items():
            if k in stem_map:
                new_samples[stem_map[k]] = v
            else:
                new_samples[k] = v
    data["samples"] = new_samples

    # 1. Check for Annotation Recomputation Bypass
    config = state.get("config", {})
    gene_name = config.get("general", {}).get("gene_name", "")
    gtfs = tuple(sorted(config.get("data_selection", {}).get("selected_gtfs", [])))
    full_range = config.get("genome_annotations", {}).get("show_full_range", False)
    
    # We only cache if 'Full Range' is UNselected, because full range discovery
    # depends on the actual reads found in BAMs, which can change.
    ann_cache_key = None
    if not full_range:
        ann_cache_key = (gene_name, gtfs, False)
        
    with state["lock"]:
        if (ann_cache_key is not None and 
            state.get("last_ann_cache_key") == ann_cache_key and 
            state.get("processed_annotations") is not None):
            log_safe(state, "Skipping redundant genome annotation recomputation (Cached).")
            return
        state["last_ann_cache_key"] = ann_cache_key

    ann_list = data.get("annotations", [])
    log_safe(state, f"Found {len(ann_list)} annotations.")
    
    if ann_list:
        df_ann = pd.DataFrame(ann_list)
        # Debug: list unique features
        unique_features = df_ann['feature'].unique()
        log_safe(state, f"DEBUG: Unique annotation features found: {unique_features}")
        
        # Exons and CDS for IGV-style UTR/coding distinction
        mask = df_ann["feature"].isin(["exon", "CDS", "five_prime_utr", "three_prime_utr", "stop_codon", "start_codon", "ncRNA", "non_coding_exon", "lnc_RNA"])
        df_raw_exons = df_ann[mask].copy()
        
        df_flat = None # Flat model for core logic
        if not df_raw_exons.empty:
            # Group by transcript and merge all overlapping/contiguous segments into true physical exons
            flat_records = []
            for tid, group in df_raw_exons.groupby('transcript_id'):
                intervals = sorted(list(zip(group['start'].astype(int), group['end'].astype(int))))
                if not intervals: continue
                
                merged = [list(intervals[0])]
                for c_start, c_end in intervals[1:]:
                    if c_start <= merged[-1][1]:
                        merged[-1][1] = max(merged[-1][1], c_end)
                    else:
                        merged.append([c_start, c_end])
                
                template = group.iloc[0]
                for m_start, m_end in merged:
                    row = template.copy()
                    row['start'], row['end'], row['feature'] = m_start, m_end, 'exon'
                    flat_records.append(row)
            
            df_flat = pd.DataFrame(flat_records).sort_values(['transcript_id', 'start'])
            
            # Recalculate physical exon numbers based on strand to fix the bug
            # for reverse strand transcripts in non-CDS areas.
            recalculated_flat = []
            for tid, group in df_flat.groupby('transcript_id', sort=False):
                strand = group['strand'].iloc[0]
                if strand == '-':
                    # For reverse strand, the first exon is at the highest coordinate
                    group = group.sort_values('start', ascending=False)
                else:
                    group = group.sort_values('start', ascending=True)
                
                group['exon_number'] = [str(i+1) for i in range(len(group))]
                recalculated_flat.append(group)
            
            if recalculated_flat:
                df_flat = pd.concat(recalculated_flat).sort_values(['transcript_id', 'start'])
            
            # Use df_flat as base for refined visualization model (splitting into UTR/CDS)
            df_exons = df_flat.copy()
            
            # Infer UTR and ncRNA status
            transcripts_with_cds = df_ann[df_ann['feature'] == 'CDS']['transcript_id'].unique()
            cds_bounds = {}
            for tid in transcripts_with_cds:
                t_cds = df_ann[(df_ann['transcript_id'] == tid) & (df_ann['feature'] == 'CDS')]
                cds_bounds[tid] = (t_cds['start'].min(), t_cds['end'].max())
            
            refined_exons = []
            for _, row in df_exons.iterrows():
                tid = row['transcript_id']
                if tid not in transcripts_with_cds:
                    row_copy = row.copy()
                    row_copy['feature'] = "ncRNA"
                    refined_exons.append(row_copy)
                    continue
                
                c_min, c_max = cds_bounds[tid]
                strand = row['strand']
                e_start, e_end = row['start'], row['end']
                
                # Part 1: Segment before CDS
                if e_start < c_min:
                    utr_end = min(e_end, c_min)
                    part = row.copy()
                    part['start'], part['end'] = e_start, utr_end
                    part['feature'] = "5'UTR" if strand == '+' else "3'UTR"
                    refined_exons.append(part)
                
                # Part 2: Segment during CDS
                coding_start = max(e_start, c_min)
                coding_end = min(e_end, c_max)
                if coding_start < coding_end:
                    part = row.copy()
                    part['start'], part['end'] = coding_start, coding_end
                    part['feature'] = "Exon (Coding)"
                    refined_exons.append(part)
                
                # Part 3: Segment after CDS
                if e_end > c_max:
                    utr_start = max(e_start, c_max)
                    part = row.copy()
                    part['start'], part['end'] = utr_start, e_end
                    part['feature'] = "3'UTR" if strand == '+' else "5'UTR"
                    refined_exons.append(part)
            
            df_exons = pd.DataFrame(refined_exons) if refined_exons else pd.DataFrame(columns=df_exons.columns)
            
        df_cds = df_ann[df_ann['feature'] == 'CDS'].copy()
        if not df_cds.empty:
            # Sort to handle strand-specific ranking
            df_plus = df_cds[df_cds['strand'] != '-'].sort_values(['transcript_id', 'start'])
            df_minus = df_cds[df_cds['strand'] == '-'].sort_values(['transcript_id', 'start'], ascending=[True, False])
            
            if not df_plus.empty:
                df_plus['cds_rank'] = (df_plus.groupby('transcript_id').cumcount() + 1).astype(str)
            if not df_minus.empty:
                df_minus['cds_rank'] = (df_minus.groupby('transcript_id').cumcount() + 1).astype(str)
            
            df_cds = pd.concat([df_plus, df_minus])
        
        if not df_exons.empty and df_flat is not None:
            known_introns = set()
            all_introns = set()
            known_sites = set()
            predicted_sites = set()
            
            # CORE LOGIC: Use Flat Model (df_flat) for introns and site assessment
            for tid, group in df_flat.groupby('transcript_id'):
                tid_str = str(tid)
                is_curated = tid_str.startswith("NM_") or tid_str.startswith("NR_")
                
                starts = group['start'].tolist()
                ends = group['end'].tolist()
                
                for i in range(len(starts) - 1):
                    # In the flat model, every gap between entries is a real intron
                    intron = (int(ends[i]), int(starts[i+1]))
                    all_introns.add(intron)
                    if is_curated:
                        known_introns.add(intron)
                        known_sites.add(int(ends[i]))
                        known_sites.add(int(starts[i+1]))
                    elif tid_str.startswith("XM_") or tid_str.startswith("XR_"):
                        predicted_sites.add(int(ends[i]))
                        predicted_sites.add(int(starts[i+1]))

            gene_col = 'gene_name' if 'gene_name' in df_exons.columns else 'gene_id' if 'gene_id' in df_exons.columns else None
            
            if gene_col:
                # IMPORTANT: Fill NaNs so they aren't dropped by groupby
                df_exons[gene_col] = df_exons[gene_col].fillna(df_exons['gene_id']).fillna("unknown")
                
                gene_data = []
                t_map = {}
                current_y = 0
                
                # Sort by start to have a consistent order
                for g_id, g_df in df_exons.sort_values('start').groupby(gene_col, sort=False):
                    transcripts = g_df['transcript_id'].unique()
                    g_min_start = g_df['start'].min()
                    g_max_end = g_df['end'].max()
                    strand = g_df['strand'].iloc[0]
                    
                    gene_y_start = current_y
                    
                    if str(g_id) == "unknown":
                        # Space efficient packing for unknown RNAs (e.state. fragments, small RNAs)
                        t_info = []
                        for tid in transcripts:
                            sub = g_df[g_df['transcript_id'] == tid]
                            t_info.append({
                                'id': tid,
                                'start': sub['start'].min(),
                                'end': sub['end'].max()
                            })
                        # Sort transcripts by start position for greedy packing
                        t_info.sort(key=lambda x: x['start'])
                        
                        rows = [] # Stores end coordinates of the last transcript in each row
                        for info in t_info:
                            placed = False
                            for row_idx, row_end in enumerate(rows):
                                # Use a small gap (e.g. 1000bp) for visual separation
                                if info['start'] > row_end + 1000:
                                    t_map[info['id']] = current_y + row_idx
                                    rows[row_idx] = info['end']
                                    placed = True
                                    break
                            if not placed:
                                t_map[info['id']] = current_y + len(rows)
                                rows.append(info['end'])
                        
                        gene_y_end = current_y + len(rows) - 1
                        current_y = current_y + len(rows)
                    else:
                        # Standard gene block: one line per transcript
                        for t_id in transcripts:
                            t_map[t_id] = current_y
                            current_y += 1
                        gene_y_end = current_y - 1
                    
                    arrow = "▶" if strand == "+" else "◀"
                    prefix = arrow * 3 + " "
                    suffix = " " + arrow * 3
                    
                    gene_data.append({
                        'text': prefix + str(g_id) + suffix,
                        'x': (g_min_start + g_max_end) / 2,
                        'y': gene_y_end + 0.6
                    })
                    current_y += 0.5
                
                df_exons['y'] = df_exons['transcript_id'].map(t_map)
                if not df_cds.empty:
                    df_cds['y'] = df_cds['transcript_id'].map(t_map)
                else:
                    df_cds['y'] = pd.Series(dtype='float64')
                df_gene_labels = pd.DataFrame(gene_data)
                if df_gene_labels.empty:
                    df_gene_labels = pd.DataFrame(columns=['text', 'x', 'y'])
            else:
                transcripts = df_exons['transcript_id'].unique()
                t_map = {t: i for i, t in enumerate(transcripts)}
                df_exons['y'] = df_exons['transcript_id'].map(lambda x: t_map.get(x, 0))
                if not df_cds.empty:
                    df_cds['y'] = df_cds['transcript_id'].map(lambda x: t_map.get(x, 0))
                else:
                    df_cds['y'] = pd.Series(dtype='float64')
                df_gene_labels = pd.DataFrame(columns=['text', 'x', 'y'])
            
            df_t = df_exons.groupby('transcript_id').agg({'start': 'min', 'end': 'max', 'y': 'first', 'strand': 'first'}).reset_index()

            # Intron markers and segments
            intron_markers = []
            intron_segments = []
            for tid, group in df_exons.groupby('transcript_id'):
                y = t_map.get(tid, 0)
                strand = group['strand'].iloc[0]
                reference = group['reference'].iloc[0] if 'reference' in group.columns else ""
                gene_id = group[gene_col].iloc[0] if gene_col else "unknown"
                starts = group['start'].tolist()
                ends = group['end'].tolist()
                for i in range(len(starts) - 1):
                    i_start = int(ends[i])
                    i_end = int(starts[i+1])
                    if i_end > i_start:
                        # Intron Segment for hover
                        intron_segments.append({
                            'start': i_start, 'end': i_end, 'y': y,
                            'transcript_id': tid, 'strand': strand,
                            'gene_id': gene_id, 'length': i_end - i_start,
                            'reference': reference
                        })
                        
                        length = i_end - i_start
                        num_markers = min(100, max(1, length // 1000))
                        step = length / (num_markers + 1)
                        for m_idx in range(1, num_markers + 1):
                            pos = i_start + m_idx * step
                            intron_markers.append({
                                'x': pos, 'y': y, 
                                'angle': -math.pi/2 if strand == '+' else math.pi/2,
                                'transcript_id': tid,
                                'strand': strand
                            })
            df_markers = pd.DataFrame(intron_markers) if intron_markers else pd.DataFrame(columns=['x', 'y', 'angle', 'transcript_id', 'strand'])
            df_introns = pd.DataFrame(intron_segments) if intron_segments else pd.DataFrame(columns=['start', 'end', 'y', 'transcript_id', 'strand', 'gene_id', 'length', 'reference'])

            with state["lock"]:
                state["known_introns"] = known_introns
                state["all_introns"] = all_introns
                state["known_sites"] = known_sites
                state["predicted_sites"] = predicted_sites
                state["processed_annotations"] = df_exons
                state["flat_annotations"] = df_flat
                state["processed_cds"] = df_cds
                state["processed_transcripts"] = df_t
                state["processed_gene_labels"] = df_gene_labels
                state["processed_markers"] = df_markers
                state["processed_introns"] = df_introns
        else:
            with state["lock"]:
                state["known_introns"] = set()
                state["all_introns"] = set()
                state["processed_annotations"] = None
                state["flat_annotations"] = None
                state["processed_cds"] = None
                state["processed_transcripts"] = None
                state["processed_gene_labels"] = None
    else:
        with state["lock"]:
            state["processed_annotations"] = None
            state["flat_annotations"] = None
            state["processed_cds"] = None
            state["processed_transcripts"] = None
            state["processed_gene_labels"] = None

def categorize_junction(state, start, end):
    """
    Categorizes a junction based on transcripts in the workspace.
    - Curated: at least one curated transcript contains this exact intron.
    - Predicted: not curated and at least one arbitrary transcript contains this exact intron.
    - Novel: Otherwise.
    """
    try:
        start_val = int(start)
        end_val = int(end)
    except (ValueError, TypeError):
        return "novel"

    known_i = state.get("known_introns", set())
    all_i = state.get("all_introns", set())
    
    j = (start_val, end_val)
    if j in known_i:
        return "curated"
    if j in all_i:
        return "predicted"
            
    return "novel"

def get_marked_sets(state, analysis_samples, mark_reqs, min_reads):
    """Calculates which junction spans should be highlighted based on comparative rules."""
    marked_sets = {"curated": set(), "predicted": set(), "novel": set()}
    for r_type in ["curated", "predicted", "novel"]:
        presence_req = mark_reqs[r_type]["presence"]
        absence_req = mark_reqs[r_type]["absence"]
        shared_set = set()
        
        if presence_req:
            per_sample_sets = []
            for name in presence_req:
                raw_name = strip_id(name)
                s = analysis_samples.get(raw_name)
                if s:
                    spans = []
                    for row_j in s.get("junction_spans", []):
                        if row_j.get("reads", 0) < min_reads: continue
                        j_start, j_end = int(row_j['start']), int(row_j['end'])
                        j_t = categorize_junction(state, j_start, j_end)
                        if j_t == r_type: spans.append((j_start, j_end))
                    per_sample_sets.append(set(spans))
            if per_sample_sets: shared_set = set.intersection(*per_sample_sets)
            
        if shared_set and absence_req:
            for name in absence_req:
                raw_name = strip_id(name)
                s = analysis_samples.get(raw_name)
                if s:
                    spans = [r for r in s.get("junction_spans", []) if r.get("reads", 0) >= min_reads]
                    j_set = set((int(r['start']), int(r['end'])) for r in spans)
                    shared_set.difference_update(j_set)
        marked_sets[r_type] = shared_set
    return marked_sets

def calculate_global_ranges(analysis_data, selected_samples_names, normalize=False, cohort_metadata=None):
    """Calculates global min/max x and max y across all selected samples."""
    global_x_min = analysis_data.get("min_x")
    global_x_max = analysis_data.get("max_x")
    global_y_max = 10.0
    
    samples = analysis_data.get("samples", {})
    
    if global_x_min is None or global_x_max is None or global_x_min >= global_x_max:
        global_x_min = None; global_x_max = None
        for name in selected_samples_names:
            raw_name = strip_id(name)
            s_data = samples.get(raw_name)
            if s_data:
                positions = s_data.get("positions", [])
                if positions:
                    lo, hi = positions[0], positions[-1]
                    global_x_min = lo if global_x_min is None else min(global_x_min, lo)
                    global_x_max = hi if global_x_max is None else max(global_x_max, hi)

    for name in selected_samples_names:
        raw_name = strip_id(name)
        s_data = samples.get(raw_name)
        if s_data:
            depths = s_data.get("depths", [])
            if depths:
                import numpy as np
                current_max = float(np.max(depths))
                
                if normalize and cohort_metadata:
                    num_samples = 1
                    if raw_name in cohort_metadata:
                        num_samples = cohort_metadata[raw_name].get("num_samples", 1)
                    current_max = current_max / max(1.0, float(num_samples))
                
                global_y_max = max(global_y_max, current_max)
                
    return global_x_min, global_x_max, global_y_max
