# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /anaconda3/lib/python3.6/site-packages/cmake/data/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /anaconda3/lib/python3.6/site-packages/cmake/data/CMake.app/Contents/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c/build

# Include any dependencies generated for this target.
include CMakeFiles/tiny-aes.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/tiny-aes.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tiny-aes.dir/flags.make

CMakeFiles/tiny-aes.dir/aes.c.o: CMakeFiles/tiny-aes.dir/flags.make
CMakeFiles/tiny-aes.dir/aes.c.o: ../aes.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/tiny-aes.dir/aes.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/tiny-aes.dir/aes.c.o   -c /Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c/aes.c

CMakeFiles/tiny-aes.dir/aes.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/tiny-aes.dir/aes.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c/aes.c > CMakeFiles/tiny-aes.dir/aes.c.i

CMakeFiles/tiny-aes.dir/aes.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/tiny-aes.dir/aes.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c/aes.c -o CMakeFiles/tiny-aes.dir/aes.c.s

# Object files for target tiny-aes
tiny__aes_OBJECTS = \
"CMakeFiles/tiny-aes.dir/aes.c.o"

# External object files for target tiny-aes
tiny__aes_EXTERNAL_OBJECTS =

libtiny-aes.a: CMakeFiles/tiny-aes.dir/aes.c.o
libtiny-aes.a: CMakeFiles/tiny-aes.dir/build.make
libtiny-aes.a: CMakeFiles/tiny-aes.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library libtiny-aes.a"
	$(CMAKE_COMMAND) -P CMakeFiles/tiny-aes.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tiny-aes.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tiny-aes.dir/build: libtiny-aes.a

.PHONY : CMakeFiles/tiny-aes.dir/build

CMakeFiles/tiny-aes.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tiny-aes.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tiny-aes.dir/clean

CMakeFiles/tiny-aes.dir/depend:
	cd /Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c /Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c /Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c/build /Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c/build /Users/ChengLuhan/Desktop/Units/FIT3143/Assignment/Assignment_2/fit3143_assignment2/code/tiny-AES-c/build/CMakeFiles/tiny-aes.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/tiny-aes.dir/depend

