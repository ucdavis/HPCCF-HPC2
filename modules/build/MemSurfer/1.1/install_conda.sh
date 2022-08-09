#!/bin/bash
set -e

source ./module.sh

module load conda3/4.13.0
source activate MemSurfer-$VERSION

pushd $SRC_DIR
python -m pip install --no-deps --force-reinstall dist/memsurfer-*-linux_x86_64.whl
popd
