#!/bin/bash
mkdir results
# init vcpkg repository
# 检查是否已定义 VCPKG_ROOT 环境变量
if [ -n "$VCPKG_ROOT" ]; then
  echo "VCPKG_ROOT is already set to: $VCPKG_ROOT"
else
  git submodule update --init --depth 1 --recursive vcpkg

  # 验证 bootstrap 脚本存在
  if [ ! -f "vcpkg/bootstrap-vcpkg.sh" ]; then
    echo "Error: bootstrap-vcpkg.sh missing in vcpkg directory!"
  fi

  # 保存初始目录
  original_dir=$(pwd)

  # 执行 bootstrap 初始化
  echo "Initializing vcpkg..."
  cd vcpkg || {
    echo "Error: Failed to enter vcpkg directory." >&2
    cd "$original_dir" # 确保回到初始目录
  }

  # 尝试执行 bootstrap 脚本
  chmod +x bootstrap-vcpkg.sh
  if ./bootstrap-vcpkg.sh; then
    # 成功时设置环境变量
    export VCPKG_ROOT="$(pwd)"
    echo -e "\nVCPKG_ROOT successfully set to: $VCPKG_ROOT"
    echo "To persist this configuration, add following to your shell profile:"
    echo "  export VCPKG_ROOT=\"$VCPKG_ROOT\""
    echo "You can now use vcpkg in this terminal session."
  else
    # 失败时回退到初始目录并报错
    echo "Error: Bootstrapping vcpkg failed. Returning to original directory." >&2
    cd "$original_dir" # 确保退出前回到初始目录
  fi

  # 无论成功与否，最终回到初始目录
  cd "$original_dir"

# 可选: 自动集成到工作流 (取消注释生效)
# ./vcpkg/vcpkg install your-packages
fi
