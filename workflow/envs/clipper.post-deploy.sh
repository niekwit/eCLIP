#!/usr/bin/env bash
# Installs the ENCODE eCLIP peak caller CLIPper (YeoLab/clipper, Python 3) into the conda environment.
# CLIPper is not available on bioconda.
# The C++ extension has a missing return statement in PyInit_peaks(), which crashes
# (segmentation fault) when compiled with a modern GCC, so this is patched before building.
set -euo pipefail

COMMIT=8bbc3db08e6bc438edb9aa11df8ed3b8f4719367

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

git clone --quiet https://github.com/YeoLab/clipper.git "$TMP/clipper"
cd "$TMP/clipper"
git checkout --quiet "$COMMIT"

sed -i 's|^  PyModule_Create(&peaksmodule);|  return PyModule_Create(\&peaksmodule);|' clipper/src/peaksmodule.cc
grep -q 'return PyModule_Create(&peaksmodule);' clipper/src/peaksmodule.cc

# Python 3.7 adds -B compiler_compat to the link command, which does not work with the conda toolchain
export LDSHARED="gcc -pthread -shared"

pip install --no-deps .
