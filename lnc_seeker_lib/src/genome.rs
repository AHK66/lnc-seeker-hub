// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::io::{self};
use std::fs::File;
use std::collections::HashMap;
use noodles::fasta;
use fasta::io::IndexedReader;
use noodles::core::Region;

pub struct GenomeProvider {
    reader: Option<IndexedReader<fasta::io::BufReader<File>>>,
    index: Option<fasta::fai::Index>,
    // Multi-reference cache: reference name -> sequence
    cache: HashMap<String, fasta::record::Sequence>,
}

impl GenomeProvider {
    pub fn new(path: &Option<String>) -> io::Result<Self> {
        let mut provider = Self {
            reader: None,
            index: None,
            cache: HashMap::new(),
        };

        if let Some(p) = path {
            let index = fasta::fai::fs::read(format!("{}.fai", p))?;
            let reader = fasta::io::indexed_reader::Builder::default()
                .build_from_path(p)?;
            provider.reader = Some(reader);
            provider.index = Some(index);
            println!("Rust Genome: Initialized indexed FASTA from {}", p);
        }

        Ok(provider)
    }

    pub fn get_sequence(&mut self, reference: &str) -> io::Result<Option<&fasta::record::Sequence>> {
        let Some(reader) = &mut self.reader else {
            return Ok(None);
        };

        if self.cache.contains_key(reference) {
            return Ok(self.cache.get(reference));
        }

        // Load new reference sequence
        println!("Rust Genome: Loading sequence for {}", reference);
        let region: Region = reference.parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let record = reader.query(&region)?;
        
        self.cache.insert(reference.to_string(), record.sequence().clone());
        
        Ok(self.cache.get(reference))
    }

    pub fn clear_cache(&mut self) {
        self.cache.clear();
    }

    /// Fetches a small slice of the sequence. 
    /// Note: This still relies on the current chromosome being cached for efficiency.
    pub fn get_range(&mut self, reference: &str, start_0: usize, end_0: usize) -> io::Result<Option<&[u8]>> {
        let Some(seq) = self.get_sequence(reference)? else {
            return Ok(None);
        };

        let seq_bytes = seq.as_ref();
        if start_0 < seq_bytes.len() {
            let actual_end = end_0.min(seq_bytes.len());
            Ok(Some(&seq_bytes[start_0..actual_end]))
        } else {
            Ok(None)
        }
    }
}
