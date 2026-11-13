cmake_minimum_required(VERSION 3.14)
project(cbeta)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

find_package(GSL REQUIRED)
LIST(APPEND EXTLIBS ${GSL_LIBRARIES})

find_package(ROOT REQUIRED COMPONENTS MathCore MathMore)
LIST(APPEND EXTLIBS ${ROOT_LIBRARIES})
include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/)
file(GLOB DIRSRC ${PROJECT_SOURCE_DIR}/*.c*)
LIST(APPEND SOURCE_FILES ${DIRSRC})

add_library(BSG SHARED ${SOURCE_FILES})
target_include_directories(BSG SYSTEM PRIVATE ${GSL_INCLUDE_DIRS})
target_compile_features(BSG PUBLIC cxx_std_14)
if(ROOT_CXX_STANDARD)
    target_compile_features(BSG PRIVATE cxx_std_${ROOT_CXX_STANDARD})
endif()
target_link_libraries(BSG ${EXTLIBS})
