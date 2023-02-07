# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.5)

# Set the project name
project (numerical_methods)

#include Eigen
include_directories(/usr/include/eigen3)



# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    src/euler_method.cpp
    src/euler_class_main.cpp
)

# Add an executable with the above sources
add_executable(euler_main src/euler_main.cpp)
add_executable(taylor_methods src/taylor_methods.cpp)


add_executable(main ${SOURCES})

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories(main
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)