cmake_minimum_required(VERSION 3.21.3 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

project(itur_p452)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10.0")
    message(FATAL_ERROR "Insufficient gcc version")
  endif()
endif()

set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)


set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++20 -O2 -g -D_GLIBCXX_DEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g") 
#set(CMAKE_CXX_FLAGS_MINSIZEREL, "-Os -DNDEBUG")
#set(CMAKE_CXX_FLAGS_RELEASE, "-O4 -DNDEBUG")

# control where the static and shared libraries are built so that on windows
# we don't need to tinker with the path to run the executable
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

if(UNIX AND NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
elseif(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWIN32_LEAN_AND_MEAN")
# elseif(MSVC)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
endif()

option(BUILD_SHARED_LIBS "build using shared libraries" ON)

find_package(Threads)
find_package(GTest REQUIRED)
find_package(GSL REQUIRED)

add_compile_definitions(CMAKE_ROOT_SRC_DIR="${CMAKE_SOURCE_DIR}/ituModels/itu_linux")
add_compile_definitions(CMAKE_CLEARAIR_SRC_DIR="${CMAKE_SOURCE_DIR}/ITUR_P452")

enable_testing()

add_subdirectory(ituModels/itu_linux/Common)
add_subdirectory(ituModels/itu_linux/GasModel)
add_subdirectory(ITUR_P452)
add_subdirectory(ClutterModel_P452_17)

