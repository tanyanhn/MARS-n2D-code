# MyToolchain.cmake
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_EXE_LINKER_FLAGS "-fuse-ld=lld")
# add_link_options("-fuse-ld=lld")

# 设置 PKG_CONFIG_PATH
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib/pkgconfig")
include_directories("${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include")
message("include_dir ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include")
# set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}:${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share")