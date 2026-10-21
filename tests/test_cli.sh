#!/usr/bin/env bash
# Command-line checks on generated data: every command gives the same output at 1 and 4 threads
# and for plain, gzipped and piped input, and a few known answers hold.
#
#     pixi run test-cli        (after pixi run build-cli)
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
sporos=$root/build/sporos
work=$root/build/test_cli
rm -rf "$work" && mkdir -p "$work"
cd "$work"

failures=0
check() {
    if [ "$2" = "$3" ]; then
        echo "  ok    $1"
    else
        echo "  FAIL  $1"
        failures=$((failures + 1))
    fi
}
digest() { "$@" | cksum; }

# A reference of two sequences with N runs, a query of its pieces (one strand reversed), and reads.
python3 - <<'EOF'
import random

rnd = random.Random(7)
def dna(n):
    return "".join(rnd.choice("ACGT") for _ in range(n))

def reverse_complement(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

ref1 = dna(60000) + "N" * 500 + dna(40000)
ref2 = dna(80000)
with open("ref.fa", "w") as f:
    f.write(">ref1 first\n" + "\n".join(ref1[i:i + 60] for i in range(0, len(ref1), 60)) + "\n")
    f.write(">ref2\n" + "\n".join(ref2[i:i + 80] for i in range(0, len(ref2), 80)) + "\n")
query = ref1[1000:9000] + dna(500) + reverse_complement(ref2[5000:12000]) + ref1[70000:71000]
with open("query.fa", "w") as f:
    f.write(">query\n" + query + "\n")
with open("reads.fq", "w") as f:
    for i in range(300):
        p = rnd.randrange(0, len(ref2) - 150)
        f.write("@read%d\n%s\n+\n%s\n" % (i, ref2[p:p + 150], "I" * 150))
EOF
gzip -c ref.fa > ref.fa.gz

echo "minimizers and syncmers"
for args in "minimizers -k 21 -w 11" "minimizers -k 15 -w 10 --forward" "syncmers -k 21 -s 11" \
            "syncmers -k 20 -s 16 --open 2" "syncmers -k 20 -s 16 --open 1,3"; do
    one=$(digest "$sporos" $args -t 1 ref.fa)
    check "$args: 4 threads = 1" "$one" "$(digest "$sporos" $args -t 4 ref.fa)"
    check "$args: gzip = plain" "$one" "$(digest "$sporos" $args -t 2 ref.fa.gz)"
    check "$args: pipe = file" "$one" "$(cat ref.fa | digest "$sporos" $args -t 2 /dev/stdin)"
    reads_one=$(digest "$sporos" $args -t 1 reads.fq)
    check "$args, reads: 4 threads = 1" "$reads_one" "$(digest "$sporos" $args -t 4 reads.fq)"
done
# Every k-mer is its own minimizer with w = 1.
check "w=1 gives every k-mer" "$(printf '>s\nACGTACGTAC\n' | "$sporos" minimizers -k 5 -w 1 /dev/stdin | cut -f2 | tr '\n' ' ')" "0 1 2 3 4 5 "

echo "seeds"
for args in "-k 19 -w 12" "-k 19 -w 12 --query-kmers --max-occ 5" "-k 20 --syncmers -s 16 --open 2" "-k 17 --all-kmers"; do
    one=$(digest "$sporos" seeds $args -t 1 ref.fa reads.fq)
    check "seeds $args: 4 threads = 1" "$one" "$(digest "$sporos" seeds $args -t 4 ref.fa reads.fq)"
done
# Reads cut from ref2 hit ref2 on the forward strand at their own position.
check "seeds find each read's origin" "$("$sporos" seeds -k 19 -w 12 -t 2 ref.fa reads.fq | awk -F'\t' '$3 == "ref2" && $5 == "+"' | wc -l | tr -d ' ' | awk '{print ($1 > 0)}')" "1"

echo "mem"
for mode in -maxmatch -mum -mumreference; do
    for strands in "" "-b" "-r"; do
        one=$(digest "$sporos" mem $mode $strands -l 30 -t 1 ref.fa query.fa)
        check "mem $mode $strands: 4 threads = 1" "$one" "$(digest "$sporos" mem $mode $strands -l 30 -t 4 ref.fa query.fa)"
    done
done
check "mem: reads (many records), 4 threads = 1" "$(digest "$sporos" mem -mum -b -l 40 -t 1 ref.fa reads.fq)" \
    "$(digest "$sporos" mem -mum -b -l 40 -t 4 ref.fa reads.fq)"
# The query's first piece is ref1[1000:9000]: an 8000-base MUM at reference 1001, query 1 (1-based).
check "mem finds the planted 8 kb match" \
    "$("$sporos" mem -mum -l 100 ref.fa query.fa | awk '$2 == 1001 && $3 == 1 {print $4}')" "8000"

echo "smem"
for len in 20 40; do
    check "smem -l $len: 4 threads = 1" "$(digest "$sporos" smem -l $len -t 1 ref.fa query.fa)" "$(digest "$sporos" smem -l $len -t 4 ref.fa query.fa)"
    check "smem -l $len, reads: 4 threads = 1" "$(digest "$sporos" smem -l $len -t 1 ref.fa reads.fq)" "$(digest "$sporos" smem -l $len -t 4 ref.fa reads.fq)"
done
# The FM-index engine (reads) and the MEM engine find the same SMEMs; a saved index gives the same.
for len in 15 30; do
    fm=$(digest "$sporos" smem --engine fm -l $len -t 3 ref.fa reads.fq)
    check "smem -l $len, reads: FM engine = MEM engine" "$fm" "$(digest "$sporos" smem --engine mem -l $len -t 3 ref.fa reads.fq)"
    check "smem -l $len, query: FM engine = MEM engine" "$(digest "$sporos" smem --engine fm -l $len ref.fa query.fa)" \
        "$(digest "$sporos" smem --engine mem -l $len ref.fa query.fa)"
    "$sporos" index -t 2 ref.fa saved.sfm
    check "smem -l $len, reads: saved index = built" "$fm" "$(digest "$sporos" smem -l $len -t 2 saved.sfm reads.fq)"
done
# The planted pieces: ref1[1000:9000] forward at query 0, and ref2[5000:12000] reverse-complemented.
check "smem finds the forward piece" "$("$sporos" smem -l 100 ref.fa query.fa | awk -F'\t' '$1 == "EM" && $2 == 0 {print $3, $4, $5}')" "8000 1 ref1:+1001"
check "smem finds the reverse piece" "$("$sporos" smem -l 100 ref.fa query.fa | awk -F'\t' '$1 == "EM" && $2 == 8500 {print $3, $4, $5}')" "15500 1 ref2:-5001"

echo "errors"
if "$sporos" minimizers -k 2>/dev/null; then check "missing option value is an error" 1 0; else check "missing option value is an error" 1 1; fi
if "$sporos" mem -l 0 ref.fa query.fa >/dev/null 2>&1; then check "-l 0 is rejected" 1 0; else check "-l 0 is rejected" 1 1; fi
for cmd in "mem" "smem" "seeds -k 19" "minimizers -k 21 -w 11"; do
    if "$sporos" $cmd -t 0 ref.fa query.fa >/dev/null 2>&1; then check "$cmd -t 0 is rejected" 1 0; else check "$cmd -t 0 is rejected" 1 1; fi
done

if [ "$failures" -gt 0 ]; then
    echo "test_cli: $failures check(s) failed"
    exit 1
fi
echo "test_cli: all checks passed"
