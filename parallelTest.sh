#!/bin/bash

tests="${1:-./build/test/TestMARSn2D}"
shift
conditions=("$@")

# 默认多条件配置
if [ ${#conditions[@]} -eq 0 ]; then
    conditions=("[Order4]" "[Graph41]")  # AND逻辑需同时满足的默认条件
fi

cond_str="${conditions[*]}"

echo "提取的测试名称："
test_list=$("$tests" --list-tests)
echo "$test_list" | awk -v cond_str="$cond_str" '
BEGIN {
    split(cond_str, conds, " ")
}
/^      \[/ {
    name = $1
    matched = 1
    for (i in conds) {
        if (conds[i] != "" && index(name, conds[i]) == 0) {
            matched = 0
            break
        }
    }
    if (matched) {
        print name
    }
}'

echo "开始运行测试..."
physical_cores=$(lscpu | awk '/^Core\(s\) per socket:/ {c=$4} /^Socket\(s\):/ {s=$2} END {print c*s}')

echo "$test_list" | awk -v cond_str="$cond_str" '
BEGIN {
    split(cond_str, conds, " ")
}
/^      \[/ {
    name = $1
    matched = 1
    for (i in conds) {
        if (conds[i] != "" && index(name, conds[i]) == 0) {
            matched = 0
            break
        }
    }
    if (matched) {
        print name
    }
}' | parallel -j "$physical_cores" "$tests {}" 1>&1