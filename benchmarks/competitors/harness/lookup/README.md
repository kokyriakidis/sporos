# Seed-index lookup harnesses

Seeding throughput of minimap2 and strobealign, measured the way `benchmarks/bench_index.mojo` measures
sporos: build an index of the reference, then on one thread seed every read and look its seeds up, timing
only the read loop (reads in memory, best of 3).

- `mm2_lookup.c`: minimap2's `mm_sketch` and `mm_idx_get` per read, copying out every hit of each
  minimizer with at most `max_occ` hits. Build against a minimap2 checkout:
  `gcc -O3 -march=native -I$M mm2_lookup.c $M/libminimap2.a -lz -lm -lpthread -o mm2_lookup`;
  run `mm2_lookup ref.fa reads.fa k w max_occ`.
- `sa_bench/`: strobealign's `randstrobes_query` and `find_hits` on both orientations, as its chainer
  calls them, with its default parameters for 150 bp reads; then again expanding every hit to its
  reference positions. `Cargo.toml` points at a strobealign checkout in `~/sporos-bench/competitors`;
  build with `RUSTFLAGS="-C target-cpu=native" cargo build --release`, run `sa_bench ref.fa reads.fa`.
