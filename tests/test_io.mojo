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


def _system(command: String) -> Int32:
    var bytes = List[UInt8](command.as_bytes())
    bytes.append(0)
    return external_call["system", Int32](bytes.unsafe_ptr())


def test_fastq_header_at_end() raises:
    var r = parse_sequences(String("@r1\nACGT\n+\nIIII\n@r2").as_bytes())
    assert_equal(len(r), 2)
    assert_true(_is(r, 0, "ACGT") and _is(r, 1, ""))
    assert_true(r.names[1] == "r2")
    assert_equal(len(parse_sequences(String("@r").as_bytes())), 1)


def test_pipe_input() raises:
    # A FIFO reports size 0: reading must still run to end of file.
    var fifo = String("build/test_io_fifo")
    assert_equal(_system("rm -f " + fifo + " && mkfifo " + fifo), 0)
    assert_equal(_system("(printf '>a\\nACGTACGT\\n>b\\nGG\\n' > " + fifo + " &)"), 0)
    var r = read_fasta(fifo)
    assert_equal(len(r), 2)
    assert_true(_is(r, 0, "ACGTACGT") and _is(r, 1, "GG"))


def test_bgzf_then_gzip_members() raises:
    if _system("command -v bgzip > /dev/null") != 0:
        return  # needs htslib's bgzip
    assert_equal(_system("printf '>a\\nACGT\\n' > build/test_io_a.fa && printf '>b\\nGGCC\\n' > build/test_io_b.fa"), 0)
    assert_equal(_system("(bgzip -c build/test_io_a.fa; gzip -c build/test_io_b.fa) > build/test_io_mixed.fa.gz"), 0)
    assert_equal(_system("(gzip -c build/test_io_a.fa; bgzip -c build/test_io_b.fa) > build/test_io_mixed2.fa.gz"), 0)
    for path in ["build/test_io_mixed.fa.gz", "build/test_io_mixed2.fa.gz"]:
        var r = read_fasta(path)
        assert_equal(len(r), 2)
        assert_true(_is(r, 0, "ACGT") and _is(r, 1, "GGCC"))


def test_parallel_parse_matches_serial() raises:
    # Over 16 MB, so read_fasta parses (and reads) on several threads.
    var path = String("build/test_io_big.fa")
    with open(path, "w") as f:
        var line = String("ACGTNacgtRYACGTACGTTTGACCAGTAGGCATACGATCGATCCGATAGCTAGCTTAGCAGT")
        for r in range(400):
            f.write(">seq", r, " description ", r * 7, "\n" if r % 3 != 0 else "\r\n")
            if r % 50 == 7:
                continue  # an empty record
            for i in range(800 + (r * 37) % 500):
                f.write(line, "\r\n" if (r + i) % 5 == 0 else "\n")
                if i % 97 == 0:
                    f.write("\n")  # a blank line
    var serial = read_fasta(path, 1)
    var parallel = read_fasta(path, 4)
    assert_true(len(serial.data) > 1 << 24)
    assert_equal(len(parallel), len(serial))
    assert_equal(len(parallel.data), len(serial.data))
    for r in range(len(serial)):
        assert_true(parallel.names[r] == serial.names[r])
        assert_equal(parallel.starts[r], serial.starts[r])
        assert_equal(parallel.lengths[r], serial.lengths[r])
    for i in range(len(serial.data)):
        if parallel.data[i] != serial.data[i]:
            raise Error(String(t"byte {i} differs"))


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
