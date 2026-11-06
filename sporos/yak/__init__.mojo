"""yak-compatible k-mer counting and hifiasm's homopolymer-compressed minimizer seeds.

hifiasm (Cheng et al. 2021) finds candidate overlaps between long reads from minimizers of the
homopolymer-compressed reads, with k up to 63 (51 by default), counted by yak's algorithm. These kernels give
the same answers as its C code, bit for bit, and are the seeding of the `hcm` package (hifiasmCandidatesMojo):

- `nt4`, `Kmer`, `yak_hash64_64` (`hash`): yak's 2-bit code, its four-plane rolling k-mers and its hash; the hash
  also eight keys at once (`yak_hash64_64_lanes`).
- `hpc_kmer_hashes`, `KmerCounts` (`count`): every compressed k-mer's hash, counted exactly or through yak's
  blocked Bloom filter (`BloomFilters`, `bloom`), which keeps the k-mers seen once out of the tables.
- `analyze_count` -> `Peaks` (`histogram`): the count histogram's homozygous and heterozygous peaks.
- `FilterTable` (`filter`): the k-mers counted at least a cut-off, as hifiasm's high-occurrence filter.
- `sketch_read` -> `Sketch` of `Minimizer` (`sketch`): hifiasm's minimizers (`mz1_ha_sketch`), ordered by the
  filter's count, then the hash, with its re-seeding of high-occurrence stretches (`select_mz_h`); `sketch_all`
  sketches every read in parallel.
- `PositionTable` (`table`): every occurrence of each minimizer kept (2 to 4094 occurrences), by hash.
"""

from .bloom import DEFAULT_BF_SHIFT, BloomFilters
from .count import KmerCounts, hpc_kmer_hashes
from .filter import FilterTable
from .hash import Kmer, nt4, yak_hash64_64, yak_hash64_64_lanes
from .histogram import Peaks, analyze_count
from .sketch import Minimizer, Sketch, mzcmp, select_mz_h, sketch_read
from .table import PositionTable, position_pos, position_rev, position_rid, position_span, sketch_all
