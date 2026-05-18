// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::collections::HashMap;

use noodles::sam;

/// Build a map of reference sequence name -> index from a SAM header.
pub fn build_name_to_index(header: &sam::Header) -> HashMap<String, usize> {
    let mut out_name_to_index: HashMap<String, usize> = HashMap::new();
    for (i, (name, _rs)) in header.reference_sequences().iter().enumerate() {
        out_name_to_index.insert(String::from_utf8_lossy(name.as_ref()).into_owned(), i);
    }
    out_name_to_index
}

/// Build a map of reference sequence index -> name from a SAM header.
pub fn build_index_to_name(header: &sam::Header) -> HashMap<usize, String> {
    let mut index_to_name: HashMap<usize, String> = HashMap::new();
    for (i, (name, _rs)) in header.reference_sequences().iter().enumerate() {
        index_to_name.insert(i, String::from_utf8_lossy(name.as_ref()).into_owned());
    }
    index_to_name
}
