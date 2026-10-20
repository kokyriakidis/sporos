#!/usr/bin/env bash
# Builds sporos's C API (see c/sporos.h) into build/c: libsporos for the platform's oldest CPU, or
# `target-cpu`, beside the Mojo runtime libraries it loads and the header.
#
#     pixi run build-c [target-cpu]
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
# The file is not c/sporos.mojo, which would be a module of the package's own name importing itself.
BUNDLE_NAME=sporos bash "$root/scripts/bundle.sh" --library "$root/c/sporos_c.mojo" "$root/build/c" "$@"
cp "$root/c/sporos.h" "$root/build/c/"
