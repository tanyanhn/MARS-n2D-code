#!/bin/bash
compiler=${1:-"gcc"}
cmake -S . -B build/ -DCMAKE_BUILD_TYPE=RELEASE --preset=$compiler