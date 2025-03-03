#!/bin/bash
source compile.sh
# make -C build test -j32
ctest --test-dir build -j32