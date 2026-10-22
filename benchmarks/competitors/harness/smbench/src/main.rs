use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec, Seq, PackedSeq};
use std::time::Instant;

fn read_fasta(path: &str) -> Vec<u8> {
    let data = std::fs::read(path).unwrap();
    let mut out = Vec::with_capacity(data.len());
    for line in data.split(|&c| c == b'\n') {
        if line.first() == Some(&b'>') { continue; }
        out.extend(line.iter().filter(|c| !c.is_ascii_whitespace()).map(|c| c.to_ascii_uppercase()));
    }
    out
}

fn run_mode(mode: &str, seq: PackedSeq, k: usize, w: usize, out: &mut Vec<u32>) {
    out.clear();
    match mode {
        "fwd" => { simd_minimizers::minimizers(k, w).run(seq, out); }
        "can" => { simd_minimizers::canonical_minimizers(k, w).run(seq, out); }
        "closed-fwd" => { simd_minimizers::closed_syncmers(k, w).run(seq, out); }
        "closed-can" => { simd_minimizers::canonical_closed_syncmers(k, w).run(seq, out); }
        "open-can" => { simd_minimizers::canonical_open_syncmers(k, w).run(seq, out); }
        "open-fwd" => { simd_minimizers::open_syncmers(k, w).run(seq, out); }
        _ => panic!(),
    }
}

fn bench(name: &str, mode: &str, pv: &PackedSeqVec, k: usize, w: usize, threads: usize, reps: usize, report: bool) {
    let seq = pv.as_slice();
    let n = seq.len();
    let mut best = f64::MAX;
    let mut count = 0usize;
    for _ in 0..reps {
        let t = Instant::now();
        if threads == 1 {
            let mut out = Vec::new();
            run_mode(mode, seq, k, w, &mut out);
            count = out.len();
        } else {
            // split into `threads*4` chunks with overlap k+w-2, process via scoped threads
            let parts = threads * 4;
            let chunk = n.div_ceil(parts);
            let next = std::sync::atomic::AtomicUsize::new(0);
            let total = std::sync::atomic::AtomicUsize::new(0);
            std::thread::scope(|s| {
                for _ in 0..threads {
                    s.spawn(|| {
                        let mut out = Vec::new();
                        loop {
                            let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            if i >= parts { break; }
                            let st = i * chunk;
                            let en = ((i + 1) * chunk + k + w - 2).min(n);
                            if st >= en || en - st < k + w - 1 { continue; }
                            run_mode(mode, seq.slice(st..en), k, w, &mut out);
                            total.fetch_add(out.len(), std::sync::atomic::Ordering::Relaxed);
                        }
                    });
                }
            });
            count = total.into_inner();
        }
        let el = t.elapsed().as_secs_f64();
        best = best.min(el);
    }
    if !report {
        return;
    }
    println!("{name:10} {mode:11} k={k:2} w={w:2} T={threads} n={n:>11} out={count:>10} density={:.4} time={:.4}s  {:.3} ns/bp  {:.3} Gbp/s",
        count as f64 / n as f64 * (w as f64 + 1.0) / 2.0 * 2.0 / (w as f64 + 1.0) * 1.0, best, best * 1e9 / n as f64, n as f64 / best / 1e9);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let threads: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(1);
    let (name, pv) = match args[1].as_str() {
        "random" => ("random100M", PackedSeqVec::random(100_000_000)),
        path => {
            let ascii = read_fasta(path);
            let t = Instant::now();
            let pv = PackedSeqVec::from_ascii(&ascii);
            let el = t.elapsed().as_secs_f64();
            eprintln!("pack {} bp ascii->2bit: {:.3}s ({:.3} ns/bp)", ascii.len(), el, el * 1e9 / ascii.len() as f64);
            (if path.contains("ecoli") { "ecoli" } else { "chm13chr1" }, pv)
        }
    };
    let reps: usize = std::env::var("REPS").ok().map(|s| s.parse().unwrap()).unwrap_or(7);
    // One untimed run first, as the sporos harnesses do: threads, pages and clocks warmed before any timing.
    bench(name, "fwd", &pv, 21, 11, threads, 1, false);
    for &(k, w) in &[(21usize, 11usize), (31, 19)] {
        for mode in ["fwd", "can", "closed-fwd", "closed-can"] {
            bench(name, mode, &pv, k, w, threads, reps, true);
        }
    }
    // strobealign default: K=20, s=16 => simd-minimizers k=16, w=5, open canonical
    bench(name, "open-fwd", &pv, 16, 5, threads, reps, true);
    bench(name, "open-can", &pv, 15, 5, threads, reps, true);
}
