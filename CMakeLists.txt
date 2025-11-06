# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and programming language
project(MyProjectName LANGUAGES CXX)

# Optionally set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add executable (change main.cpp to your source file)
add_executable(MyExecutable main.cpp)