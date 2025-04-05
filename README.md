```markdown
# MARS: Multi YinSet Simulation in 2D Space
<!-- [![Build Status](https://img.shields.io/badge/build-passing-brightgreen)]()
[![License](https://img.shields.io/badge/license-MIT-blue)]() -->
Front Tracking method for multi phases Interface Tracking problems in 2D space.

## 🚀 System Requirements

### Essential Environment
- Ubuntu 22.04 LTS (other Linux distributions may require dependency adjustments)
- C++20 compatible compilers:
  - GCC 11+ or
  - Clang 14+

### Dependency Installation
```bash
# Base build tools
sudo apt update && sudo apt install -y \
    cmake ninja-build pkg-config \
    git curl autoconf libtool

# Development libraries
sudo apt install -y \
    liblapacke-dev \
    libomp-dev  # For Clang users: libomp-dev
```

## 🛠️ Build Guide

### 1. Environment Initialization
```bash
source init.sh  # Sets environment variables and paths
```

### 2. Configuration & Compilation
```bash
# Select compiler (gcc/clang)
source releaseconfigure.sh gcc  # or clang

# Parallel compilation
source compile.sh
```

## 🧪 Test Verification
```bash
# Run full test suite
source parallelTest.sh
```