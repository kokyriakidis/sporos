# sporos

Fast DNA seeding in Mojo: (w, k)-minimizers, open/closed syncmers, MEMs and MUMs. The goal is to be the
fastest library for each, by a wide margin; rivals are simd-minimizers, minimap2, MUMmer4, copMEM2, mumemto,
BWA (SMEMs) and GraphAligner's seeders.

## Mojo skills

Modular's Mojo skills live in `.claude/skills` and correct pretrained Mojo, which is out of date. Hooks in
`.claude/settings.json` inject them: the full `mojo-syntax` skill and `references/idiomatic-mojo.md` at session
start (also after `/clear` and compaction) and into every subagent, and a reminder on every prompt
(`SPOROS_SKILLS_EVERY_PROMPT=full` injects the full text on every prompt instead). Follow them over pretrained
habits:

- `mojo-syntax` before writing, editing or reviewing any Mojo, with `references/idiomatic-mojo.md`.
- `mojo-gpu-fundamentals` before touching GPU code; `mojo-python-interop` before Mojo-Python interop.
- Include the skills rule in every subagent prompt.

`pixi run sync-skills` refreshes them from upstream; `.claude/skills/UPSTREAM` records the commit.

## Working here

- The public API is what `sporos/__init__.mojo` exports (the README's Library section shows it); keep the
  two in step when adding or renaming anything public.
- Build and test: `pixi run test` (every kernel against the brute-force oracles in `tests/reference.mojo`);
  format: `pixi run format` (120 columns).
- Toolchain: one Mojo nightly, pinned exactly in `pixi.toml` (1.2.0.dev2026102005). Move it on purpose: change the
  pin, `pixi lock`, build and test on xaros, then re-run the benchmarks sporos-only.
- Pointers handed out by owning types must carry the owner's origin (see `CodeBuffer.unsafe_ptr`); an
  untracked pointer lets the owner be destroyed early. Buffers that tasks reach through untracked pointers
  are explicit `alloc`/`dealloc` pairs: `_ = list^` after the tasks does not keep a `List` alive. Prefer tasks
  generic over the origins of the data they read (`cli/mem_command.mojo`'s `_FormatTask`): tracked pointers keep
  the owners alive while the task runs.
- `unsafe_strided_load` on `UInt8` returns wrong lanes on the current nightly; gather explicitly.
- In a `TestSuite` file, `String(StringSlice(unsafe_from_utf8=span))` fails to compile (recursive
  `CStringSpan.__init__`); compare bytes instead. `std.subprocess` trips the same bug.
- Benchmark with min-of-N timings; other jobs share the machine.
- Finished work merges straight into `main` and is pushed; no pull requests.
- The paper's text lives in `tools/paper/` (`main_text.py`, `supplement.py`, `refs.py`); its numbers come from
  `benchmarks/competitors/results/` through `tables.py`. `pixi run paper` rewrites `paper/`, an Overleaf project,
  from them, and `pixi run paper-sync` pushes it (the `overleaf` remote; the token stays in the keychain); with
  `pixi run paper-hook` installed, every commit that changes `paper/` pushes it. Edit
  the sources, never the generated `.tex`; type no number into the text that `tables.py` can compute.
