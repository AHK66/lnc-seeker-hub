// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use serde::{Serialize, Deserialize};

/// Genomic read name edit operations for differential encoding.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum EditOp {
    Match(u8),         // Number of bytes to copy from predecessor
    Insert(Vec<u8>),   // New bytes to insert
    Delete(u8),        // Bytes to skip in predecessor
    Substitute(Vec<u8>), // Replacement bytes (advances both)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CompressionMode {
    None,
    Huffman,
    Zstd,
}

/// A block of compressed header data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompressedHeaders {
    pub mode: CompressionMode,
    pub data: Vec<u8>,
    pub symbol_lengths: Option<Vec<u8>>, // Only for Huffman
    pub num_names: usize,
    pub use_substitutes: bool,
}

// Alphabet:
// 0-255: Literal byte
// 256-511: Match(len) where len = sym - 256
// 512-767: Delete(len) where len = sym - 512
// 768-1023: Substitute(len) where len = sym - 768 (Optional)
// 768 or 1024: EndOfName
const SYMBOL_SUBSTITUTE_START: u16 = 768;
const MAX_SYMBOLS_STD: usize = 1025;
const MAX_SYMBOLS_COMPACT: usize = 769;
const LUT_BITS: u8 = 11;

pub struct BitWriter {
    buffer: Vec<u8>,
    bit_buf: u64,
    bits_count: u8,
}

impl BitWriter {
    pub fn new() -> Self {
        Self { buffer: Vec::with_capacity(65536), bit_buf: 0, bits_count: 0 }
    }

    #[inline(always)]
    pub fn write_bits(&mut self, code: u32, len: u8) {
        if len == 0 { return; }
        // Ensure space in buffer
        if self.bits_count + len > 64 {
            while self.bits_count >= 8 {
                self.buffer.push((self.bit_buf >> 56) as u8);
                self.bit_buf <<= 8;
                self.bits_count -= 8;
            }
        }
        self.bit_buf |= (code as u64) << (64 - self.bits_count - len);
        self.bits_count += len;
    }

    pub fn finish(mut self) -> Vec<u8> {
        while self.bits_count > 0 {
            self.buffer.push((self.bit_buf >> 56) as u8);
            self.bit_buf <<= 8;
            if self.bits_count >= 8 {
                self.bits_count -= 8;
            } else {
                self.bits_count = 0;
            }
        }
        self.buffer
    }
}

pub struct FastBitReader<'a> {
    data: &'a [u8],
    pos: usize,
    bit_buf: u64,
    bits_count: u8,
}

impl<'a> FastBitReader<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        let mut r = Self { data, pos: 0, bit_buf: 0, bits_count: 0 };
        r.refill();
        r
    }

    #[inline(always)]
    fn refill(&mut self) {
        while self.bits_count <= 56 && self.pos < self.data.len() {
            self.bit_buf |= (self.data[self.pos] as u64) << (56 - self.bits_count);
            self.pos += 1;
            self.bits_count += 8;
        }
    }

    #[inline(always)]
    pub fn peek_bits(&mut self, n: u8) -> u32 {
        (self.bit_buf >> (64 - n)) as u32
    }

    #[inline(always)]
    pub fn consume_bits(&mut self, n: u8) {
        self.bit_buf <<= n;
        self.bits_count -= n;
        if self.bits_count < 32 {
            self.refill();
        }
    }
}

pub fn compute_name_diff(old: &[u8], new: &[u8], use_substitutes: bool) -> Vec<EditOp> {
    let mut ops = Vec::new();
    let n = old.len();
    let m = new.len();
    let mut i = 0; 
    let mut j = 0; 

    while i < n && j < m {
        let mut match_len = 0;
        while i + match_len < n && j + match_len < m && old[i + match_len] == new[j + match_len] {
            match_len += 1;
        }
        if match_len > 0 {
            let mut rem = match_len;
            while rem > 0 {
                let take = rem.min(255) as u8;
                ops.push(EditOp::Match(take));
                rem -= take as usize;
            }
            i += match_len;
            j += match_len;
            if i >= n && j >= m { break; }
        }

        let mut found_sub = false;
        if use_substitutes {
            for sub_len in 1..=16 {
                if i + sub_len < n && j + sub_len < m {
                    let mut re_align = 0;
                    while i + sub_len + re_align < n && j + sub_len + re_align < m 
                        && old[i + sub_len + re_align] == new[j + sub_len + re_align] {
                        re_align += 1;
                    }
                    if re_align >= 4 { 
                        ops.push(EditOp::Substitute(new[j..j + sub_len].to_vec()));
                        i += sub_len;
                        j += sub_len;
                        found_sub = true;
                        break;
                    }
                }
            }
        }

        if !found_sub {
            let old_rem = n - i;
            if old_rem > 0 {
                let mut rem = old_rem;
                while rem > 0 {
                    let take = rem.min(255) as u8;
                    ops.push(EditOp::Delete(take));
                    rem -= take as usize;
                }
                i = n;
            }
            if j < m {
                ops.push(EditOp::Insert(new[j..].to_vec()));
                j = m;
            }
        }
    }
    
    if i < n {
        let mut rem = n - i;
        while rem > 0 {
            let take = rem.min(255) as u8;
            ops.push(EditOp::Delete(take));
            rem -= take as usize;
        }
    }
    if j < m {
        ops.push(EditOp::Insert(new[j..].to_vec()));
    }
    ops
}

pub fn names_to_symbols(names: &[String], use_substitutes: bool) -> Vec<u16> {
    let mut symbols = Vec::with_capacity(names.len() * 10);
    let mut prev: Vec<u8> = Vec::new();
    let end_symbol = if use_substitutes { 1024 } else { 768 };

    for name in names {
        let current = name.as_bytes();
        let ops = compute_name_diff(&prev, current, use_substitutes);
        for op in ops {
            match op {
                EditOp::Match(len) => symbols.push(256 + len as u16),
                EditOp::Delete(len) => symbols.push(512 + len as u16),
                EditOp::Substitute(data) => {
                    symbols.push(SYMBOL_SUBSTITUTE_START + data.len() as u16);
                    for &b in &data { symbols.push(b as u16); }
                }
                EditOp::Insert(data) => {
                    for &b in &data { symbols.push(b as u16); }
                }
            }
        }
        symbols.push(end_symbol);
        prev = current.to_vec();
    }
    symbols
}

pub fn symbols_to_names(symbols: &[u16], num_names: usize, use_substitutes: bool) -> Vec<String> {
    let mut names = Vec::with_capacity(num_names);
    let mut prev = Vec::new();
    let mut current_name = Vec::with_capacity(128);
    let mut old_pos = 0;
    let end_symbol = if use_substitutes { 1024 } else { 768 };

    for &sym in symbols {
        if sym == end_symbol {
            names.push(unsafe { String::from_utf8_unchecked(current_name.clone()) });
            std::mem::swap(&mut prev, &mut current_name);
            current_name.clear();
            old_pos = 0;
        } else if sym < 256 {
            current_name.push(sym as u8);
        } else if sym < 512 {
            let mlen = (sym - 256) as usize;
            if old_pos + mlen <= prev.len() {
                current_name.extend_from_slice(&prev[old_pos..old_pos + mlen]);
                old_pos += mlen;
            }
        } else if sym < 768 {
            old_pos += (sym - 512) as usize;
        } else if use_substitutes && sym < 1024 {
            old_pos += (sym - 768) as usize;
        }
    }
    names
}

pub fn compress_header_block(names: &[String], mode: CompressionMode, use_substitutes: bool) -> CompressedHeaders {
    let symbols = names_to_symbols(names, use_substitutes);
    
    match mode {
        CompressionMode::None => {
            let mut data = Vec::with_capacity(symbols.len() * 2);
            for &s in &symbols { data.extend_from_slice(&s.to_le_bytes()); }
            CompressedHeaders { mode, data, symbol_lengths: None, num_names: names.len(), use_substitutes }
        },
        CompressionMode::Zstd => {
            let mut raw = Vec::with_capacity(symbols.len() * 2);
            for &s in &symbols { raw.extend_from_slice(&s.to_le_bytes()); }
            let data = zstd::encode_all(&raw[..], 3).expect("Zstd compression failed");
            CompressedHeaders { mode, data, symbol_lengths: None, num_names: names.len(), use_substitutes }
        },
        CompressionMode::Huffman => {
            let num_symbols = if use_substitutes { MAX_SYMBOLS_STD } else { MAX_SYMBOLS_COMPACT };
            let mut freqs = vec![0usize; num_symbols];
            for &s in &symbols { freqs[s as usize] += 1; }
            let lengths = generate_canonical_lengths(&freqs);
            let mut codes = vec![(0u32, 0u8); num_symbols];
            let mut max_len = 0;
            for &l in &lengths { if l > max_len { max_len = l; } }

            let mut next_code = 0u32;
            for len in 1..=max_len {
                for (sym, &l) in lengths.iter().enumerate() {
                    if l == len {
                        codes[sym] = (next_code, len);
                        next_code += 1;
                    }
                }
                next_code <<= 1;
            }

            let mut writer = BitWriter::new();
            for &sym in &symbols {
                let (c, l) = codes[sym as usize];
                writer.write_bits(c, l);
            }
            CompressedHeaders {
                mode,
                data: writer.finish(),
                symbol_lengths: Some(lengths),
                num_names: names.len(),
                use_substitutes,
            }
        }
    }
}

pub fn decompress_header_block(compressed: &CompressedHeaders) -> Vec<String> {
    if compressed.num_names == 0 { return Vec::new(); }

    match compressed.mode {
        CompressionMode::None => {
            let symbols: Vec<u16> = compressed.data.chunks_exact(2)
                .map(|c| u16::from_le_bytes([c[0], c[1]])).collect();
            symbols_to_names(&symbols, compressed.num_names, compressed.use_substitutes)
        },
        CompressionMode::Zstd => {
            let decompressed = zstd::decode_all(&compressed.data[..]).expect("Zstd decompression failed");
            let symbols: Vec<u16> = decompressed.chunks_exact(2)
                .map(|c| u16::from_le_bytes([c[0], c[1]])).collect();
            symbols_to_names(&symbols, compressed.num_names, compressed.use_substitutes)
        },
        CompressionMode::Huffman => {
            let symbol_lengths = compressed.symbol_lengths.as_ref().expect("Huffman requires symbol lengths");
            let mut reader = FastBitReader::new(&compressed.data);
            let mut max_len = 0;
            for &l in symbol_lengths { if l > max_len { max_len = l; } }
            
            let mut lut = [0u32; 1 << LUT_BITS];
            let mut long_codes = Vec::new();
            let mut next_code = 0u32;
            for len in 1..=max_len {
                for (sym, &l) in symbol_lengths.iter().enumerate() {
                    if l == len {
                        if len <= LUT_BITS {
                            let shift = LUT_BITS - len;
                            let start = (next_code << shift) as usize;
                            let end = start + (1 << shift);
                            let val = ((sym as u32) << 6) | (len as u32);
                            for i in start..end { lut[i] = val; }
                        } else {
                            long_codes.push((len, next_code, sym as u16));
                        }
                        next_code += 1;
                    }
                }
                next_code <<= 1;
            }

            let mut names = Vec::with_capacity(compressed.num_names);
            let mut prev = Vec::new();
            let mut current_name = Vec::with_capacity(128);
            let mut old_pos = 0;
            let mut names_decoded = 0;
            let end_symbol = if compressed.use_substitutes { 1024 } else { 768 };

            while names_decoded < compressed.num_names {
                let peek = reader.peek_bits(LUT_BITS);
                let entry = lut[peek as usize];
                let sym = if entry != 0 {
                    let s = (entry >> 6) as u16;
                    reader.consume_bits((entry & 0x3F) as u8);
                    s
                } else {
                    let mut res = None;
                    for &(l, c, s) in &long_codes {
                        if reader.peek_bits(l) == c {
                            reader.consume_bits(l);
                            res = Some(s);
                            break;
                        }
                    }
                    if let Some(s) = res { s } else { break; }
                };

                if sym == end_symbol {
                    names.push(unsafe { String::from_utf8_unchecked(current_name.clone()) });
                    std::mem::swap(&mut prev, &mut current_name);
                    current_name.clear();
                    old_pos = 0;
                    names_decoded += 1;
                } else if sym < 256 {
                    current_name.push(sym as u8);
                } else if sym < 512 {
                    let mlen = (sym - 256) as usize;
                    if old_pos + mlen <= prev.len() {
                        current_name.extend_from_slice(&prev[old_pos..old_pos + mlen]);
                        old_pos += mlen;
                    }
                } else if sym < 768 {
                    old_pos += (sym - 512) as usize;
                } else if compressed.use_substitutes && sym < 1024 {
                    old_pos += (sym - 768) as usize;
                }
            }
            names
        }
    }
}

pub fn generate_canonical_lengths(freqs: &[usize]) -> Vec<u8> {
    use std::collections::BinaryHeap;
    use std::cmp::Ordering;

    #[derive(Eq, PartialEq)]
    struct Node { freq: usize, symbol: Option<u16>, left: Option<Box<Node>>, right: Option<Box<Node>> }
    impl Ord for Node { fn cmp(&self, other: &Self) -> Ordering { other.freq.cmp(&self.freq) } }
    impl PartialOrd for Node { fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) } }

    let num_symbols = freqs.len();
    let mut heap = BinaryHeap::new();
    for (i, &f) in freqs.iter().enumerate() {
        if f > 0 { heap.push(Node { freq: f, symbol: Some(i as u16), left: None, right: None }); }
    }
    if heap.is_empty() { return vec![0; num_symbols]; }
    if heap.len() == 1 {
        let node = heap.pop().unwrap();
        let mut lengths = vec![0; num_symbols];
        lengths[node.symbol.unwrap() as usize] = 1;
        return lengths;
    }
    while heap.len() > 1 {
        let left = heap.pop().unwrap();
        let right = heap.pop().unwrap();
        heap.push(Node { freq: left.freq + right.freq, symbol: None, left: Some(Box::new(left)), right: Some(Box::new(right)) });
    }
    let root = heap.pop().unwrap();
    let mut lengths = vec![0; num_symbols];
    fn walk(node: &Node, depth: u8, lengths: &mut [u8]) {
        if let Some(sym) = node.symbol { lengths[sym as usize] = depth; }
        else {
            if let Some(l) = &node.left { walk(l, depth + 1, lengths); }
            if let Some(r) = &node.right { walk(r, depth + 1, lengths); }
        }
    }
    walk(&root, 0, &mut lengths);
    lengths
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_compression_modes() {
        let names = vec![
            "K00188:452:HLYLWBBXX:1:1101:1446:1152".to_string(),
            "K00188:452:HLYLWBBXX:1:1101:1446:1153".to_string(),
            "K00188:452:HLYLWBBXX:1:1101:1532:1152".to_string(),
            "K00188:452:HLYLWBBXX:1:1101:1532:1154".to_string(),
            "K00188:452:HLYLWBBXX:1:1112:1532:1154".to_string(),
        ];

        let modes = vec![CompressionMode::None, CompressionMode::Huffman, CompressionMode::Zstd];

        for &mode in &modes {
            for &use_subs in &[true, false] {
                println!("Testing mode: {:?}, substitutes: {}", mode, use_subs);
                let compressed = compress_header_block(&names, mode, use_subs);
                assert_eq!(compressed.num_names, names.len());
                assert_eq!(compressed.mode, mode);
                assert_eq!(compressed.use_substitutes, use_subs);

                let decompressed = decompress_header_block(&compressed);
                assert_eq!(decompressed, names);
            }
        }
    }

    #[test]
    fn test_empty_block() {
        let names: Vec<String> = Vec::new();
        let compressed = compress_header_block(&names, CompressionMode::Huffman, true);
        assert_eq!(compressed.num_names, 0);
        let decompressed = decompress_header_block(&compressed);
        assert!(decompressed.is_empty());
    }
}

