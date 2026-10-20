#!/usr/bin/env bash
# Builds the Python package (see python/sporos) into build/python: the package's Python, its
# extension `_sporos` for the platform's oldest CPU or `target-cpu`, and the Mojo runtime libraries it
# loads, so `PYTHONPATH=build/python python3 -c "import sporos"` works without pixi or Mojo.
#
#     pixi run build-python [target-cpu]
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
out=$root/build/python/sporos
mkdir -p "$out"
bash "$root/scripts/bundle.sh" --library "$root/python/_sporos.mojo" "$out" "$@"
# Python imports an extension by its module's name, the bundle names a library `lib<name>`.
for library in "$out"/lib_sporos.so "$out"/lib_sporos.dylib; do
    [ -e "$library" ] && mv "$library" "$out/_sporos.so"
done
cp "$root/python/sporos/"*.py "$out/"
echo "$out: import sporos with $(dirname "$out") on PYTHONPATH"
