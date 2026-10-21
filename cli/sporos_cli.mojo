"""The `sporos` command line: `minimizers`, `syncmers`, `seeds`, `mem`,
`smem` and `index` (see `USAGE` in `common.mojo`, one module per command).

`sporos mem` takes MUMmer's `mummer` flags and prints its output format
(`-mumreference` is the default, as in MUMmer; non-ACGT bases never match);
`sporos smem` prints `bwa fastmap`'s.
"""

from std.sys import argv

from cli.common import USAGE
from cli.kmers_command import run_seeds
from cli.mem_command import run_mem
from cli.seeds_command import run_seed_index
from cli.smem_command import run_index, run_smem


def main() raises:
    var raw = argv()
    if len(raw) < 2:
        print(USAGE)
        return
    var command = String(raw[1])
    var args = List[String]()
    for i in range(2, len(raw)):
        args.append(String(raw[i]))
    if command == "mem":
        run_mem(args)
    elif command == "seeds":
        run_seed_index(args)
    elif command == "smem":
        run_smem(args)
    elif command == "index":
        run_index(args)
    elif command == "minimizers" or command == "syncmers":
        run_seeds(command, args)
    else:
        print(USAGE)
