#!/bin/bash
tests=${1:-"./build/test/TestMARSn2D"}
echo "提取的测试名称："
$tests --list-tests | awk '/^      \[/ {print $1}'
echo "开始运行测试..."
$tests --list-tests | awk '/^      \[/ {print $1}' | parallel -j4 "$tests "{}"" 1>&1
