#!/usr/bin/env bash
# Minimizer/syncmer benchmark: sporos vs simd-minimizers vs minimap2's mm_sketch.
#
#     bash benchmarks/competitors/run_seeds.sh
#
# Kernel timings (no FASTA parsing, no output) on 100 Mbp random DNA and CHM13 chr1:
#   sporos           benchmarks/bench_seeds.mojo (random) and bench_seeds_fasta.mojo (chr1),
#                    built with SPOROS_LANES=8 and 16; input already 2-bit coded, 1 byte/base
#   simd-minimizers  harness/smbench; input already packed 2 bits/base
#   minimap2         harness/mm2 (mm_sketch on ASCII; k=31 runs as k=28, its maximum)
# Single thread: best of 5-7 runs (thread CPU time for sporos, wall time for the others).
# 10 threads: best wall time; simd-minimizers and mm_sketch split the sequence into chunks.
# Also the end-to-end CLI: `sporos minimizers -k 21 -w 11 -t T chr1.fa > /dev/null`.
# Writes results/seeds.tsv (parsed), results/seeds_raw.txt and results/seeds_cli.tsv.
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"

CHR1=$D/chm13_chr1.fa
SMB=$BENCH/mybench/smbench/target/release/smbench
MMB=$BENCH/mybench/mm2/mmbench
TSV=$RES/seeds.tsv
RAW=$RES/seeds_raw.txt
CLI=$RES/seeds_cli.tsv
printf "tool\tinput\tmode\tk\tw\tthreads\tn\tout\tns_per_bp\tbusy_cpus\tload1\n" > "$TSV"
: > "$RAW"
{ echo "== run_seeds.sh"; machine_info; } >> "$RES/machine.txt"

# run <tool label> <input label> cmd...: raw lines to $RAW, key=value fields to $TSV.
run() {
  local tool=$1 input=$2; shift 2
  wait_idle
  local busy ld; busy=$(busy_cpus); ld=$(load1)
  echo "## $tool $input: $*" >> "$RAW"
  "$@" 2>&1 | tee -a "$RAW" | awk -v tool="$tool" -v input="$input" -v busy="$busy" -v ld="$ld" '
    /ns\/bp/ && /k=/ && /out=/ {
      line = $0; gsub(/= +/, "=", line); n = split(line, a, " ")
      mode = ""; k = ""; w = ""; t = 1; len = ""; out = ""; ns = ""
      for (i = 1; i <= n; i++) {
        if (a[i] ~ /^k=/) k = substr(a[i], 3)
        else if (a[i] ~ /^w=/) w = substr(a[i], 3)
        else if (a[i] ~ /^T=/) t = substr(a[i], 3)
        else if (a[i] ~ /^n=/) len = substr(a[i], 3)
        else if (a[i] ~ /^out=/) out = substr(a[i], 5)
        else if (a[i] ~ /^ns\/bp=/) ns = substr(a[i], 7)
        else if (a[i] == "ns/bp") ns = a[i - 1]
        else if (a[i] ~ /^(fwd|can|closed-fwd|closed-can|open-fwd|open-can)$/) mode = a[i]
      }
      if (mode == "" && $1 == "mm_sketch") mode = "can"
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.4f\t%s\t%s\n", tool, input, mode, k, w, t, len, out, ns, busy, ld
    }' | tee -a "$TSV"
}

# Single thread, 100 Mbp random DNA.
for lanes in 8 16; do
  run "sporos (lanes=$lanes)" random100M "$SBIN/bench_seeds_$lanes" 100000000
done
run simd-minimizers random100M "$SMB" random 1
run "minimap2 mm_sketch" random100M "$MMB" random 1

# Single thread, CHM13 chr1.
for lanes in 8 16; do
  run "sporos (lanes=$lanes)" chm13_chr1 "$SBIN/bench_seeds_fasta_$lanes" "$CHR1" 1
done
run simd-minimizers chm13_chr1 "$SMB" "$CHR1" 1
run "minimap2 mm_sketch" chm13_chr1 "$MMB" "$CHR1" 1

# 10 threads, CHM13 chr1 (wall time).
run "sporos (lanes=16)" chm13_chr1 "$SBIN/bench_seeds_fasta_16" "$CHR1" 10
run simd-minimizers chm13_chr1 "$SMB" "$CHR1" 10
run "minimap2 mm_sketch" chm13_chr1 "$MMB" "$CHR1" 10

# End-to-end CLI (read FASTA, encode, compute, print name\tpos lines).
printf "pair\tL\ttool\tmode\tthreads\treps\twall_s\tuser_s\tsys_s\tmax_rss_mb\tbusy_cpus\tload1\n" > "$CLI"
for T in 1 10; do
  measure "$CLI" chm13_chr1 - "sporos minimizers -k 21 -w 11" canonical "$T" -- - \
    "$SBIN/sporos" minimizers -k 21 -w 11 -t "$T" "$CHR1"
done
echo "run_seeds.sh done"
