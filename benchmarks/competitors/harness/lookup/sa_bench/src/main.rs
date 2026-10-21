//! strobealign seed lookups: index a reference with strobealign's defaults for
//! 150 bp reads, then for every read generate its randstrobes (both
//! orientations) and look them up with find_hits, as strobealign's chainer
//! does; then again expanding every hit to its reference positions.
//!   sa_bench <ref.fa> <reads.fa>
use std::io::{BufRead, BufReader};
use std::time::Instant;
use strobealign::hit::find_hits;
use strobealign::indexer::make_index;
use strobealign::io::fasta::read_ref;
use strobealign::mcsstrategy::McsStrategy;
use strobealign::seeding::{SeedingParameters, randstrobes_query};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let refseq = read_ref(&args[1]).unwrap();
    let params = SeedingParameters::new(150);
    let bits = params.syncmer.pick_bits(&refseq);
    for threads in [1, 8] {
        let t = Instant::now();
        let _ = make_index(&refseq, params.clone(), bits, 0.0002, threads);
        println!("build threads={threads}: {:.3} s", t.elapsed().as_secs_f64());
    }
    let (index, _) = make_index(&refseq, params.clone(), bits, 0.0002, 1);
    let mut reads: Vec<Vec<u8>> = vec![];
    for line in BufReader::new(std::fs::File::open(&args[2]).unwrap()).lines() {
        let line = line.unwrap();
        if !line.starts_with('>') {
            reads.push(line.into_bytes());
        }
    }
    let bases: usize = reads.iter().map(|r| r.len()).sum();
    for expand in [false, true] {
        let mut best = f64::MAX;
        let mut total = 0usize;
        for _ in 0..3 {
            total = 0;
            let mut positions: Vec<(u32, usize)> = Vec::with_capacity(1 << 16);
            let t = Instant::now();
            for read in &reads {
                let strobes = randstrobes_query(read, &index.parameters);
                for orientation in 0..2 {
                    let (hits, _) = find_hits(&strobes[orientation], &index, McsStrategy::Always, index.filter_cutoff(), 100);
                    if !expand {
                        total += hits.len();
                        continue;
                    }
                    positions.clear();
                    for hit in &hits {
                        if hit.is_filtered {
                            continue;
                        }
                        let entry = index.entry(hit.position);
                        let count = if hit.is_partial { entry.get_count_partial() } else { entry.get_count_full_forward() };
                        for i in 0..count {
                            positions.push((hit.query_start as u32, index.entry(hit.position + i).ref_start()));
                        }
                    }
                    total += positions.len();
                }
            }
            best = best.min(t.elapsed().as_secs_f64());
        }
        let what = if expand { "find_hits + every position" } else { "find_hits (one per seed)" };
        println!("lookup {what}: {:.1} Mbp/s, {:.2} per read ({best:.3} s)", bases as f64 / best / 1e6, total as f64 / reads.len() as f64);
    }
}
