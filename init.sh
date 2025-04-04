#!/bin/bash
mkdir results
# init vcpkg repository
# 检查是否已定义 VCPKG_ROOT 环境变量
if [ -n "$VCPKG_ROOT" ]; then
    echo "VCPKG_ROOT is already set to: $VCPKG_ROOT"
else
# 当任何命令返回非零时退出脚本
# set -e

# 如果 vcpkg 目录不存在
if [ ! -d "vcpkg" ]; then
    echo "VCPKG_ROOT not set and vcpkg directory not found. Setting up vcpkg..."

    # 检测是否在 Git 仓库中
    if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        echo "Detected Git repository. Initializing vcpkg submodule..."
        
        # 尝试更新现有子模块（浅克隆）
        if ! git submodule update --init --depth 1 --recursive vcpkg 2>/dev/null; then
            echo "Adding vcpkg as new submodule..."
            # 添加子模块时指定深度为1
            git submodule add --depth 1 https://github.com/Microsoft/vcpkg.git vcpkg
            # 初始化子模块（浅克隆）
            git submodule update --init --depth 1 --recursive vcpkg
        fi
    else  
        echo "Should work in a git repo..."
        # echo "Cloning vcpkg repository directly with shallow clone..."
        # # 直接克隆时指定深度为1
        # git clone --depth 1 https://github.com/Microsoft/vcpkg.git vcpkg
    fi
fi

# 验证 bootstrap 脚本存在
if [ ! -f "vcpkg/bootstrap-vcpkg.sh" ]; then
    echo "Error: bootstrap-vcpkg.sh missing in vcpkg directory!"
    exit 1
fi

# 执行 bootstrap 初始化
echo "Initializing vcpkg..."
(
    cd vcpkg
    chmod +x bootstrap-vcpkg.sh
    ./bootstrap-vcpkg.sh
)

# 设置环境变量并输出提示
export VCPKG_ROOT="$(pwd)/vcpkg"
echo -e "\nVCPKG_ROOT successfully set to: $VCPKG_ROOT"
echo "To persist this configuration, add following to your shell profile:"
echo "  export VCPKG_ROOT=\"$VCPKG_ROOT\""
echo "You can now use vcpkg in this terminal session."

# 可选: 自动集成到工作流 (取消注释生效)
# ./vcpkg/vcpkg install your-packages
fi
