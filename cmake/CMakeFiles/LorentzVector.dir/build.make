# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/Veritas_Ordo/Physics

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/Veritas_Ordo/Physics/cmake

# Include any dependencies generated for this target.
include CMakeFiles/LorentzVector.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/LorentzVector.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/LorentzVector.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/LorentzVector.dir/flags.make

CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.o: CMakeFiles/LorentzVector.dir/flags.make
CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.o: ../package/LorentzVector/LorentzVector.cpp
CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.o: CMakeFiles/LorentzVector.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/Veritas_Ordo/Physics/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.o -MF CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.o.d -o CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.o -c /home/Veritas_Ordo/Physics/package/LorentzVector/LorentzVector.cpp

CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/Veritas_Ordo/Physics/package/LorentzVector/LorentzVector.cpp > CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.i

CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/Veritas_Ordo/Physics/package/LorentzVector/LorentzVector.cpp -o CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.s

# Object files for target LorentzVector
LorentzVector_OBJECTS = \
"CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.o"

# External object files for target LorentzVector
LorentzVector_EXTERNAL_OBJECTS =

libLorentzVector.a: CMakeFiles/LorentzVector.dir/package/LorentzVector/LorentzVector.cpp.o
libLorentzVector.a: CMakeFiles/LorentzVector.dir/build.make
libLorentzVector.a: CMakeFiles/LorentzVector.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/Veritas_Ordo/Physics/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libLorentzVector.a"
	$(CMAKE_COMMAND) -P CMakeFiles/LorentzVector.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LorentzVector.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/LorentzVector.dir/build: libLorentzVector.a
.PHONY : CMakeFiles/LorentzVector.dir/build

CMakeFiles/LorentzVector.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/LorentzVector.dir/cmake_clean.cmake
.PHONY : CMakeFiles/LorentzVector.dir/clean

CMakeFiles/LorentzVector.dir/depend:
	cd /home/Veritas_Ordo/Physics/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/Veritas_Ordo/Physics /home/Veritas_Ordo/Physics /home/Veritas_Ordo/Physics/cmake /home/Veritas_Ordo/Physics/cmake /home/Veritas_Ordo/Physics/cmake/CMakeFiles/LorentzVector.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/LorentzVector.dir/depend
