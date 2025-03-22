#!/bin/bash
tests=${1:-"./build/test/TestMARSn2D"}
echo "提取的测试名称："
$tests --list-tests | awk '/^      \[/ {print $1}'
echo "开始运行测试..."
# 计算物理核心数
physical_cores=$(lscpu | awk '/^Core\(s\) per socket:/ {c=$4} /^Socket\(s\):/ {s=$2} END {print c*s}')
$tests --list-tests | awk '/^      \[/ {print $1}' | parallel -j $physical_cores "$tests "{}"" 1>&1
