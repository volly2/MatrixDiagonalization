cmake_minimum_required(VERSION 3.5)
project(MatterEffect)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)

include_directories(~/ClionProjects/MatterEffect)

add_executable(MatterEffect ${SOURCE_FILES})


