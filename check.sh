#!/bin/bash
set -e 
source init.sh

compilers=("clang" "gcc")
for compiler in "${compilers[@]}"; do
  source clear.sh
  source releaseconfigure.sh $compiler
  source test.sh

  source clear.sh
  source debugconfigure.sh $compiler
  source test.sh
  source clear.sh
done
