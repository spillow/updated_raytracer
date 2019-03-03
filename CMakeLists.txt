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

set(GLUT_ROOT_PATH libs/FreeGLUT/freeglut/freeglut)
set(GLUI_ROOT_PATH libs/glui)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(FREEGLUT_BUILD_DEMOS OFF)

add_subdirectory(${GLUT_ROOT_PATH})
add_subdirectory(${GLUI_ROOT_PATH})

target_include_directories(glui_obj PUBLIC
    ${CMAKE_SOURCE_DIR}/${GLUI_ROOT_PATH}/include)
target_include_directories(glui PUBLIC
    ${CMAKE_SOURCE_DIR}/${GLUI_ROOT_PATH}/include)
target_include_directories(glui_static PUBLIC
    ${CMAKE_SOURCE_DIR}/${GLUI_ROOT_PATH}/include)

target_include_directories(RayTracer PUBLIC
    include
    ${GLUT_ROOT_PATH}/include
    ${GLUI_ROOT_PATH}/include)
