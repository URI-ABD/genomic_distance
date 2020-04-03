extern crate bio;
extern crate rayon;
extern crate num;
// extern crate distances;

// use bstr::ext_slice::ByteSlice;
use rayon::prelude::*;
use std::io;
// use std::String;
use bio::io::fasta;

fn main() {
    let reader = fasta::Reader::new(io::stdin());
    let mut seqs: Vec<String> = Vec::new();
    for result in reader.records() {
    	let record = result.unwrap();
    	let seq = String::from_utf8_lossy(record.seq());
    	seqs.push(seq.to_string());
    	// let row = reader.records().iter().map(|r1, r2| levenshtein(r1.seq(), r2.seq()));
    }
    for s1 in &seqs {
    	let row: Vec<usize>  = seqs.par_iter().map(|s2| levenshtein(&s1, &s2)).collect();
    	print!("{:?}\n", row);
    }
}

pub fn levenshtein(a: &str, b: &str) -> usize {
    let len_a = a.chars().count();
    let len_b = b.chars().count();
    if len_a < len_b{
        return levenshtein(b, a)
    }
    // handle special case of 0 length
    if len_a == 0 {
        return len_b
    } else if len_b == 0 {
        return len_a
    }

    let len_b = len_b + 1;

    let mut pre;
    let mut tmp;

    // initialize DP table for string b
    let mut cur: Vec<usize> = (0..len_b).collect();

    // calculate edit distance
    for (i,ca) in a.chars().enumerate() {
        // get first column for this row
        pre = cur[0];
        cur[0] = i + 1;
        for (j, cb) in b.chars().enumerate() {
            tmp = cur[j + 1];
            cur[j + 1] = std::cmp::min(
                // deletion
                tmp + 1, std::cmp::min(
                // insertion
                cur[j] + 1,
                // match or substitution
                pre + if ca == cb { 0 } else { 1 }));
            pre = tmp;
        }
    }
    cur[len_b - 1]
}