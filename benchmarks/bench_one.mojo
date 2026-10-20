"""One seed kernel run for profiling: bench_one <mode> [length] (mode: fwd, can, closed-fwd, closed-can)."""

from std.sys import argv

from sporos.alphabet import CodeBuffer
from sporos.minimizers import minimizer_positions_from_codes
from sporos.syncmers import syncmer_positions_from_codes
from sporos.window import PADDING
from reference import random_dna


def main() raises:
    var args = argv()
    var mode = String(args[1])
    var n = Int(args[2]) if len(args) > 2 else 100_000_000
    var seq = random_dna(n, 5)
    var codes = CodeBuffer.encode(Span(seq), PADDING)
    var total = 0
    for _ in range(5):
        if mode == "fwd" or mode == "can":
            total += len(minimizer_positions_from_codes(codes.unsafe_ptr(), n, 21, 11, mode == "can", 1))
        else:
            total += len(syncmer_positions_from_codes(codes.unsafe_ptr(), n, 31, 21, mode == "closed-can", None, 1))
    print(mode, total)
