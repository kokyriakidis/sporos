# sporos

Fast DNA seeding in [Mojo](https://mojolang.org): (w, k)-minimizers, open and closed syncmers,
maximal exact matches (MEMs) and maximal unique matches (MUMs).

Status: early development. Every kernel is tested against a brute-force oracle; benchmarks against
simd-minimizers, minimap2, MUMmer4, copMEM2 and mumemto are in progress.

## Design

- **Minimizers and syncmers** share one engine (`sporos/window.mojo`): the sequence is split into 8
  SIMD lanes, each rolls a 32-bit ntHash (canonical = forward + reverse complement) and keeps a
  sliding-window minimum with the prefix/suffix-block method. Seed lookups are single byte shuffles
  (NEON `tbl` / x86 `pshufb`), and keys are single u32 values (hash bits above a block offset), which
  keeps leftmost tie-breaking exact.
- **MEMs** (`sporos/mems.mojo`) use sampled k-mer anchors: reference K-mers every `k1` positions,
  query K-mers every `k2`, with `k1`, `k2` coprime and `k1 * k2 <= L - K + 1`. By the Chinese remainder
  theorem each MEM is found from exactly one anchor (the one whose left extension is shorter than
  `k1 * k2`), so no deduplication pass is needed. Extensions compare 8 bytes at a time.
- **MUMs** are the MEMs not contained, in reference or query coordinates, by another MEM of the same
  query sequence: exactly the matches whose string occurs once in each.
- Threads come from a dependency-free pthreads `parallel_for` (`sporos/parallel.mojo`).

## Use

```mojo
from sporos.minimizers import minimizer_positions
from sporos.syncmers import syncmer_positions

var seq = "ACGTTGCA...".as_bytes()
var mins = minimizer_positions(seq, k=21, w=11)                  # canonical by default
var syncs = syncmer_positions(seq, k=21, s=11)                   # closed syncmers
var opens = syncmer_positions(seq, k=21, s=11, open_offset=5)    # open syncmers
```

## Develop

Requires [pixi](https://pixi.sh). `pixi run test` builds and runs the test suites; `pixi run format`
formats the sources. Coding agents: see `CLAUDE.md` (the Modular Mojo skills are vendored in
`.claude/skills` and injected by hooks).
