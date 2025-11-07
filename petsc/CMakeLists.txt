# see example at https://gitlab.com/petsc/petsc/-/blob/5fa70555f2cfa5f8527759fb2fd8b5523acdf153/share/petsc/CMakeLists.txt
cmake_minimum_required(VERSION 3.15...3.31)
set(PROJECT_NAME TsunamiPropagation)
set(LIB_NAME tsunami)
set(PACKAGE_VERSION 0.1.0)

project(
  ${PROJECT_NAME}
  VERSION ${PACKAGE_VERSION}
  LANGUAGES C)

find_package(PkgConfig REQUIRED)

pkg_search_module(PETSC REQUIRED IMPORTED_TARGET)
if(PETSC_VERSION)
    message(STATUS "PETSC Version: ${PETSC_VERSION}")
endif()

add_executable(${LIB_NAME} main.c)
target_link_libraries(${LIB_NAME} PkgConfig::PETSC)
