from std.testing import assert_equal, assert_true, TestSuite
from std.ffi import external_call

from sporos.fasta import FastaRecords, parse_sequences, read_fasta


def _is(records: FastaRecords, i: Int, expected: StringSlice) -> Bool:
    """Whether record i's sequence is `expected` (compared byte by byte: building a
    String from the span here trips a compiler bug in TestSuite files)."""
    var seq = records.seq(i)
    var want = expected.as_bytes()
    if len(seq) != len(want):
        return False
    for j in range(len(seq)):
        if seq[j] != want[j]:
            return False
    return True


def test_fasta_multiline_crlf_blank_lines() raises:
    var text = String(">a desc\r\nACGT\r\nTT\r\n\r\n>b\nGG\nCC\n>c\n")
    var r = parse_sequences(text.as_bytes())
    assert_equal(len(r), 3)
    assert_true(r.names[0] == "a")
    assert_true(_is(r, 0, "ACGTTT"))
    assert_true(_is(r, 1, "GGCC"))
    assert_true(r.names[2] == "c")
    assert_true(_is(r, 2, ""))


def test_fastq() raises:
    var text = String("@r1 x\nACGTN\n+\nIIIII\n@r2\r\nGG\r\n+r2\r\n##\r\n")
    var r = parse_sequences(text.as_bytes())
    assert_equal(len(r), 2)
    assert_true(r.names[0] == "r1")
    assert_true(_is(r, 0, "ACGTN"))
    assert_true(r.names[1] == "r2")
    assert_true(_is(r, 1, "GG"))


def test_gzip_round_trip() raises:
    var path = String("build/test_io_reads.fq")
    with open(path, "w") as f:
        f.write("@q\nACGTACGT\n+\nIIIIIIII\n")
    var command = List[UInt8]((String("gzip -kf ") + path).as_bytes())
    command.append(0)
    assert_equal(external_call["system", Int32](command.unsafe_ptr()), 0)
    var plain = read_fasta(path)
    var packed = read_fasta(path + ".gz")
    assert_equal(len(packed), 1)
    assert_true(_is(packed, 0, "ACGTACGT") and _is(plain, 0, "ACGTACGT"))


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
