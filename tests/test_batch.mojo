from std.testing import assert_equal, TestSuite

from sporos.batch import minimizers_batch, syncmers_batch
from sporos.minimizers import minimizer_positions
from sporos.syncmers import syncmer_positions
from reference import random_dna


def _bounds(n: Int) -> List[Tuple[Int, Int]]:
    """Read extents in one buffer: varied lengths (0 to ~600, some below k or
    w), with gaps."""
    var out = List[Tuple[Int, Int]]()
    var at = 0
    for i in range(n):
        var length = (i * 37) % 600 + (i % 7)
        if i % 11 == 0:
            length = i % 25
        if i % 97 == 50:
            length = 9000 + i  # long enough to be scanned alone
        out.append((at, at + length))
        at += length + 13
    return out^


def _check_same(got: Span[UInt32, _], expected: List[UInt32], label: String) raises:
    if len(got) != len(expected):
        raise Error(label + String(t": {len(got)} seeds, expected {len(expected)}"))
    for i in range(len(got)):
        if got[i] != expected[i]:
            raise Error(label + String(t": seed {i} differs"))


def test_minimizers_batch_matches_single_reads() raises:
    var all = random_dna(200_000, 3, 0.003)
    var bounds = _bounds(533)
    var spans = List[Span[UInt8, origin_of(all)]]()
    for b in bounds:
        spans.append(Span(all)[b[0] : b[1]])
    for cfg in [(21, 11), (15, 10), (31, 19), (11, 5)]:
        for canonical in [True, False]:
            for threads in [1, 3]:
                var batch = minimizers_batch(spans, cfg[0], cfg[1], canonical, threads)
                assert_equal(len(batch), len(spans))
                for i in range(len(spans)):
                    var expected = minimizer_positions(spans[i], cfg[0], cfg[1], canonical)
                    _check_same(batch.seeds(i), expected, String(t"minimizers read {i} k={cfg[0]} w={cfg[1]}"))


def test_syncmers_batch_matches_single_reads() raises:
    var all = random_dna(200_000, 9, 0.003)
    var bounds = _bounds(301)
    var spans = List[Span[UInt8, origin_of(all)]]()
    for b in bounds:
        spans.append(Span(all)[b[0] : b[1]])
    for cfg in [(21, 11), (20, 16), (15, 5)]:
        for canonical in [True, False]:
            var closed = syncmers_batch(spans, cfg[0], cfg[1], canonical)
            var t = (cfg[0] - cfg[1]) // 2
            var opened = syncmers_batch(spans, cfg[0], cfg[1], canonical, t, 2)
            for i in range(len(spans)):
                _check_same(closed.seeds(i), syncmer_positions(spans[i], cfg[0], cfg[1], canonical), "closed")
                _check_same(opened.seeds(i), syncmer_positions(spans[i], cfg[0], cfg[1], canonical, t), "open")


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
