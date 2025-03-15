#!/bin/bash
compiler=${1:-"clang"}
cmake -S . -B build/ -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_CXX_FLAGS=-DOPTNONE --preset=$compiler