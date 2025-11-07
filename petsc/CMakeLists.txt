cmake_minimum_required(VERSION 3.15...3.31)
set(PROJECT_NAME TsunamiPropagation)
set(LIB_NAME tsunami)
set(PACKAGE_VERSION 0.1.0)

project(
  ${PROJECT_NAME}
  VERSION ${PACKAGE_VERSION}
  LANGUAGES C)

add_executable(${LIB_NAME} main.c)
#target_link_libraries(${LIB_NAME} PETSC)
