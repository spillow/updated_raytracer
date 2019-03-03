cmake_minimum_required(VERSION 3.0)

project(RayTracer)

set(SRC_FILES
    src/GUI.cpp
    src/main.cpp
    src/mathobjects.cpp
    src/objects.cpp
    src/RayTracer.cpp)

set(HDR_FILES
    include/includes.h
    include/mathobjects.h)

add_executable(RayTracer ${SRC_FILES} ${HDR_FILES})

target_include_directories(RayTracer PUBLIC
    include)
