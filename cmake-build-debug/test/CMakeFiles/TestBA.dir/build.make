# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ubuntu/codebase/YinSets2D-master

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ubuntu/codebase/YinSets2D-master/cmake-build-debug

# Include any dependencies generated for this target.
include test/CMakeFiles/TestBA.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/TestBA.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/TestBA.dir/flags.make

test/CMakeFiles/TestBA.dir/TestBA.cpp.o: test/CMakeFiles/TestBA.dir/flags.make
test/CMakeFiles/TestBA.dir/TestBA.cpp.o: ../test/TestBA.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/codebase/YinSets2D-master/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/TestBA.dir/TestBA.cpp.o"
	cd /home/ubuntu/codebase/YinSets2D-master/cmake-build-debug/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestBA.dir/TestBA.cpp.o -c /home/ubuntu/codebase/YinSets2D-master/test/TestBA.cpp

test/CMakeFiles/TestBA.dir/TestBA.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestBA.dir/TestBA.cpp.i"
	cd /home/ubuntu/codebase/YinSets2D-master/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/codebase/YinSets2D-master/test/TestBA.cpp > CMakeFiles/TestBA.dir/TestBA.cpp.i

test/CMakeFiles/TestBA.dir/TestBA.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestBA.dir/TestBA.cpp.s"
	cd /home/ubuntu/codebase/YinSets2D-master/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/codebase/YinSets2D-master/test/TestBA.cpp -o CMakeFiles/TestBA.dir/TestBA.cpp.s

# Object files for target TestBA
TestBA_OBJECTS = \
"CMakeFiles/TestBA.dir/TestBA.cpp.o"

# External object files for target TestBA
TestBA_EXTERNAL_OBJECTS =

test/TestBA: test/CMakeFiles/TestBA.dir/TestBA.cpp.o
test/TestBA: test/CMakeFiles/TestBA.dir/build.make
test/TestBA: src/libbays.a
test/TestBA: test/CMakeFiles/TestBA.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/codebase/YinSets2D-master/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable TestBA"
	cd /home/ubuntu/codebase/YinSets2D-master/cmake-build-debug/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestBA.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/TestBA.dir/build: test/TestBA

.PHONY : test/CMakeFiles/TestBA.dir/build

test/CMakeFiles/TestBA.dir/clean:
	cd /home/ubuntu/codebase/YinSets2D-master/cmake-build-debug/test && $(CMAKE_COMMAND) -P CMakeFiles/TestBA.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/TestBA.dir/clean

test/CMakeFiles/TestBA.dir/depend:
	cd /home/ubuntu/codebase/YinSets2D-master/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/codebase/YinSets2D-master /home/ubuntu/codebase/YinSets2D-master/test /home/ubuntu/codebase/YinSets2D-master/cmake-build-debug /home/ubuntu/codebase/YinSets2D-master/cmake-build-debug/test /home/ubuntu/codebase/YinSets2D-master/cmake-build-debug/test/CMakeFiles/TestBA.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/TestBA.dir/depend

