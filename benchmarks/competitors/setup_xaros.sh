#!/usr/bin/env bash
# Idempotent setup of the competitor benchmark on a Linux x86-64 box (xaros).
#
# Run on the benchmark machine from the sporos checkout (mirrored there with
# `rsync -az --exclude .pixi --exclude build ./ xaros:sporos/`):
#
#     bash benchmarks/competitors/setup_xaros.sh
#
# Everything lands under $BENCH (default ~/sporos-bench):
#   env/          pixi env with bioconda mummer4, mumemto, bwa, ropebwt3 (+ e-mem if solvable)
#   competitors/  copMEM2, minimap2, simd-minimizers, minibwa, strobealign sources at pinned commits
#   mybench/      simd-minimizers (Rust) and minimap2 mm_sketch (C) harnesses, and the seed-lookup
#                 harnesses of minimap2 and strobealign (harness/lookup)
#   data/         genome pairs (see fetch_data below)
#   reads/        simulated reads from hg38 chr21 for the SMEM benchmark (sim_reads.py)
#   versions.tsv  tool versions/commits
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
BENCH=${BENCH:-$HOME/sporos-bench}
PIXI=${PIXI:-$(command -v pixi || echo "$HOME/.pixi/bin/pixi")}
mkdir -p "$BENCH"/{env,competitors,mybench,data}

# --- bioconda tools ---------------------------------------------------------
if [ ! -f "$BENCH/env/pixi.toml" ]; then
  cat > "$BENCH/env/pixi.toml" <<'EOF'
[workspace]
name = "sporos-competitors"
channels = ["conda-forge", "bioconda"]
platforms = ["linux-64"]

[dependencies]
mummer4 = "*"
mumemto = "*"
bwa = "*"
ropebwt3 = "*"
EOF
fi
(cd "$BENCH/env" && "$PIXI" install)
# E-MEM is optional: keep the env usable if bioconda has no linux build.
if ! grep -q '^e-mem' "$BENCH/env/pixi.toml"; then
  (cd "$BENCH/env" && "$PIXI" add e-mem) || echo "e-mem: not installable from bioconda, skipped"
fi
ENVBIN="$BENCH/env/.pixi/envs/default/bin"

# --- sources at pinned commits ----------------------------------------------
clone() {  # url dir commit
  local url=$1 dir=$BENCH/competitors/$2 commit=$3
  [ -d "$dir/.git" ] || git clone -q "$url" "$dir"
  if [ "$(git -C "$dir" rev-parse HEAD)" != "$commit" ]; then
    git -C "$dir" fetch -q origin "$commit" 2>/dev/null || git -C "$dir" fetch -q origin
    git -C "$dir" checkout -q "$commit"
  fi
}
clone https://github.com/wbieniec/copmem2.git copmem2 8040165d635db113ee5a71ddd279f21cf57cebc4
clone https://github.com/lh3/minimap2.git minimap2 3c28777e7e2dcc90f825de1b9f17a89cca7d4452
clone https://github.com/rust-seq/simd-minimizers.git simd-minimizers a173cbe9c2bca81f622ce9296996f694c997f071
clone https://github.com/lh3/minibwa.git minibwa 74bdd485f2517fdd322faa1fd8eb36ecf297743d
clone https://github.com/ksahlin/strobealign.git strobealign e45e291a41e850d857fef9efaea35a31406c61ee

# minibwa: its own makefile (gcc -O3).
[ -x "$BENCH/competitors/minibwa/minibwa" ] || make -C "$BENCH/competitors/minibwa" -j8 >/dev/null

# copMEM2: its makefile already uses -O3 -march=native -funroll-loops; gcc >= 13 needs <cstdint> forced in.
if [ ! -x "$BENCH/competitors/copmem2/copmem2" ]; then
  make -C "$BENCH/competitors/copmem2" CPP="g++ -include cstdint" copmem2 >/dev/null
fi

# minimap2 mm_sketch harness.
mkdir -p "$BENCH/mybench/mm2"
cp "$HERE/harness/mm2/mmbench.c" "$BENCH/mybench/mm2/mmbench.c"
MM2=$BENCH/competitors/minimap2
if [ ! -f "$MM2/libminimap2.a" ]; then
  make -C "$MM2" CFLAGS="-O3 -march=native -Wall" libminimap2.a >/dev/null
fi
if [ ! -x "$BENCH/mybench/mm2/mmbench" ] || [ "$BENCH/mybench/mm2/mmbench.c" -nt "$BENCH/mybench/mm2/mmbench" ]; then
  gcc -O3 -march=native -I"$MM2" "$BENCH/mybench/mm2/mmbench.c" "$MM2/libminimap2.a" \
    -lz -lpthread -lm -o "$BENCH/mybench/mm2/mmbench"
fi

# simd-minimizers harness (depends on ../../competitors/simd-minimizers).
mkdir -p "$BENCH/mybench/smbench/src"
cp "$HERE/harness/smbench/Cargo.toml" "$HERE/harness/smbench/Cargo.lock" "$BENCH/mybench/smbench/"
cp "$HERE/harness/smbench/src/main.rs" "$BENCH/mybench/smbench/src/main.rs"
(cd "$BENCH/mybench/smbench" && RUSTFLAGS="-C target-cpu=native" cargo build -q --release)

# Seed-lookup harnesses: minimap2's mm_sketch + mm_idx_get, strobealign's find_hits (its Cargo.toml finds the
# strobealign checkout above by a path relative to this checkout, so both must sit in the home directory).
mkdir -p "$BENCH/mybench/lookup"
gcc -O3 -march=native -I"$MM2" "$HERE/harness/lookup/mm2_lookup.c" "$MM2/libminimap2.a" -lz -lm -lpthread \
  -o "$BENCH/mybench/lookup/mm2_lookup"
(cd "$HERE/harness/lookup/sa_bench" && RUSTFLAGS="-C target-cpu=native" cargo build -q --release \
  --target-dir "$BENCH/mybench/lookup/sa_bench-target")
cp "$BENCH/mybench/lookup/sa_bench-target/release/sa_bench" "$BENCH/mybench/lookup/sa_bench"

# --- sporos -------------------------------------------------------------------
(cd "$REPO" && "$PIXI" run -q build-cli && "$PIXI" run -q build-bench-seeds && "$PIXI" run -q build-bench-index)
(cd "$REPO" && for lanes in 8 16; do
  "$PIXI" run -q mojo build -D ASSERT=none -D SPOROS_LANES=$lanes -I . -I tests -I benchmarks \
    benchmarks/bench_seeds.mojo -o build/bench_seeds_$lanes
  "$PIXI" run -q mojo build -D ASSERT=none -D SPOROS_LANES=$lanes -I . -I benchmarks \
    benchmarks/competitors/bench_seeds_fasta.mojo -o build/bench_seeds_fasta_$lanes
done)
# Snapshot the binaries so later rebuilds of the checkout don't change a run midway.
mkdir -p "$BENCH/bin"
(cd "$REPO/build" && cp sporos bench_index bench_seeds_8 bench_seeds_16 bench_seeds_fasta_8 bench_seeds_fasta_16 "$BENCH/bin/")

# --- data ---------------------------------------------------------------------
# E. coli K-12/O157 and hg38/panTro chr21 come from the sporos dev box (rsync to
# $BENCH/data) or from NCBI/UCSC:
#   ecoli_k12.fa   NC_000913.3   ecoli_o157.fa NC_002695.2
#   hg38_chr21.fa  UCSC hg38 chr21   panTro_chr21_PTRv2.fa NC_036900.1
# The scaling pair is CHM13v2.0 chr1 vs hg38 chr1.
D=$BENCH/data
fetch_ncbi() {  # accession out
  [ -s "$2" ] || curl -sfL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$1&rettype=fasta&retmode=text" -o "$2"
}
fetch_ncbi NC_000913.3 "$D/ecoli_k12.fa"
fetch_ncbi NC_002695.2 "$D/ecoli_o157.fa"
fetch_ncbi NC_036900.1 "$D/panTro_chr21_PTRv2.fa"
[ -s "$D/hg38_chr21.fa" ] || curl -sfL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz | gunzip > "$D/hg38_chr21.fa"
[ -s "$D/hg38_chr1.fa" ] || curl -sfL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz | gunzip > "$D/hg38_chr1.fa"
if [ ! -s "$D/chm13_chr1.fa" ]; then
  if [ -s "$D/chm13.fa.gz" ]; then
    gunzip -c "$D/chm13.fa.gz" | awk '/^>/{p=($1==">chr1")} p' > "$D/chm13_chr1.fa"
  else
    curl -sfL https://hgdownload.soe.ucsc.edu/goldenPath/hs1/chromosomes/chr1.fa.gz | gunzip > "$D/chm13_chr1.fa"
  fi
fi

# Reads for the SMEM benchmark (run_smems.sh), drawn from hg38 chr21 with fixed seeds.
mkdir -p "$BENCH/reads"
[ -s "$BENCH/reads/chr21_150bp_1M.fa" ] || python3 -I "$HERE/sim_reads.py" "$D/hg38_chr21.fa" 1000000 150 0.002 1 "$BENCH/reads/chr21_150bp_1M.fa"
[ -s "$BENCH/reads/chr21_10kb_4k.fa" ] || python3 -I "$HERE/sim_reads.py" "$D/hg38_chr21.fa" 4000 10000 0.002 2 "$BENCH/reads/chr21_10kb_4k.fa"

# --- versions -------------------------------------------------------------------
{
  printf "tool\tversion\n"
  printf "sporos\t%s\n" "$(git -C "$REPO" describe --always --dirty 2>/dev/null || echo mirror) (mojo $(cd "$REPO" && "$PIXI" run -q mojo --version 2>/dev/null | awk '{print $2}'))"
  printf "mummer\t%s\n" "$("$ENVBIN/mummer" --version 2>&1 | head -1)"
  printf "mumemto\t%s\n" "$(cd "$BENCH/env" && "$PIXI" list mumemto 2>/dev/null | awk '$1=="mumemto"{print $2}')"
  printf "bwa\t%s\n" "$("$ENVBIN/bwa" 2>&1 | awk '/^Version/{print $2}')"
  printf "ropebwt3\t%s\n" "$("$ENVBIN/ropebwt3" version 2>&1 | head -1)"
  [ -x "$ENVBIN/e-mem" ] && printf "e-mem\t%s\n" "$(cd "$BENCH/env" && "$PIXI" list e-mem 2>/dev/null | awk '$1=="e-mem"{print $2}')"
  printf "strobealign\tgit %s (rustc, target-cpu=native)\n" "$(git -C "$BENCH/competitors/strobealign" rev-parse --short HEAD)"
  printf "minibwa\tgit %s (gcc -O3)\n" "$(git -C "$BENCH/competitors/minibwa" rev-parse --short HEAD)"
  printf "copMEM2\tgit %s (g++ -O3 -march=native)\n" "$(git -C "$BENCH/competitors/copmem2" rev-parse --short HEAD)"
  printf "minimap2 mm_sketch\tgit %s (gcc -O3 -march=native)\n" "$(git -C "$MM2" rev-parse --short HEAD)"
  printf "simd-minimizers\tgit %s (rustc %s, target-cpu=native)\n" \
    "$(git -C "$BENCH/competitors/simd-minimizers" rev-parse --short HEAD)" "$(rustc --version | awk '{print $2}')"
  printf "gcc\t%s\n" "$(gcc -dumpfullversion)"
} > "$BENCH/versions.tsv"
cat "$BENCH/versions.tsv"
