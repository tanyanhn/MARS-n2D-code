# MyToolchain.cmake
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_EXE_LINKER_FLAGS "-fuse-ld=bfd")

# 设置 PKG_CONFIG_PATH
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib/pkgconfig")
include_directories("${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include")