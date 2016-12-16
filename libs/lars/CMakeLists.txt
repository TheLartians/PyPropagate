cmake_minimum_required(VERSION 3.1)

# ---- Macros -----

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

# ---- Project ----
project(Lars)
set(LARS_LIBRARY TRUE CACHE INTERNAL "")

# ---- Configuration variables ----
set(LARS_BUILD_EXAMPLES FALSE CACHE BOOL "build the examples included in the lars libraries")
set(LARS_BUILD_TESTS FALSE CACHE BOOL "build the tests included in the lars libraries")

# ---- Init lists ----

FILE(GLOB lars_headers "${CMAKE_CURRENT_SOURCE_DIR}/include/lars/*.h")
list(APPEND lars_include_dirs "${CMAKE_CURRENT_SOURCE_DIR}/include")

# ---- Add modules ----

SUBDIRLIST(modules ${CMAKE_CURRENT_SOURCE_DIR}/modules)

FOREACH(module ${modules})
  ADD_SUBDIRECTORY(${module})
ENDFOREACH()

add_library(Lars EXCLUDE_FROM_ALL ${lars_headers})
SET_TARGET_PROPERTIES(Lars PROPERTIES LINKER_LANGUAGE CXX)

set(Lars_INCLUDE_DIRS ${lars_include_dirs} CACHE INTERNAL "")

