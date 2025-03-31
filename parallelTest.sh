#!/bin/bash

tests=${1:-"./build/test/TestMARSn2D"}
condition=${2:-"[Graph41]"}

echo "提取的测试名称："
test_list=$("$tests" --list-tests)
echo "$test_list" | awk -v cond="$condition" '
    /^      \[/ {
        name = $1
        if (cond == "" || index(name, cond) > 0) {
            print name
        }
    }'

echo "开始运行测试..."
physical_cores=$(lscpu | awk '/^Core\(s\) per socket:/ {c=$4} /^Socket\(s\):/ {s=$2} END {print c*s}')

echo "$test_list" | awk -v cond="$condition" '
    /^      \[/ {
        name = $1
        if (cond == "" || index(name, cond) > 0) {
            print name
        }
    }' | parallel -j "$physical_cores" "$tests {}" 1>&1