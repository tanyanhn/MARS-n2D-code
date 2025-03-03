#!/bin/bash
compiler=${1:-"clang"}
cmake -S . -B build/ -DCMAKE_BUILD_TYPE=RELEASE --preset=$compiler