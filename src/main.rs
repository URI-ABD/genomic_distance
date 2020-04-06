extern crate bio;
extern crate rayon;
extern crate num;
// extern crate distances;

// use bstr::ext_slice::ByteSlice;
use rayon::prelude::*;
use std::io;
// use std::String;
use bio::io::fasta;
use bio::alignment::pairwise::banded::*;
use bio::scores::blosum62;

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
    	let row: Vec<f64>  = seqs.par_iter().map(|s2| nwm(&s1, &s2)).collect();
    	print!("{:?}\n", row);
    }
}

pub fn nwm(a: &str, b: &str) -> f64 {
	/*
	Needleman-Wunsch metric, as described by Fisher

	*/
	let gapo = -5;
	let gape = -1;
	let x = a.as_bytes();
	let y = b.as_bytes();
	let score = &blosum62;
	let k=4;
	let w=12;
	let mut aligner_x = Aligner::with_capacity(x.len(), x.len(), gapo, gape, &score,k,w);
	let mut aligner_y = Aligner::with_capacity(y.len(), y.len(), gapo, gape, &score,k,w);
	let mut aligner = Aligner::with_capacity(x.len(), y.len(), gapo, gape, &score,k,w);
	let x_y = aligner.global(x, y).score;
	let x_x = aligner_x.global(x, x).score;
	let y_y = aligner_y.global(y, y).score;
	((x_x + y_y - 2*x_y) as f64).sqrt()

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