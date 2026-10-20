"""Minimizers of many 1 kb reads, for profiling per-call cost."""

from sporos.minimizers import minimizer_positions
from reference import random_dna


def main() raises:
    var seq = random_dna(10_000_000, 11)
    var count = 0
    for _ in range(3):
        for r in range(10_000):
            count += len(minimizer_positions(Span(seq)[r * 1000 : (r + 1) * 1000], 21, 11, True, 1))
    print(count)
