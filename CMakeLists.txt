
cmake_minimum_required(VERSION 2.8.3)

GET_FILENAME_COMPONENT(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
PROJECT (${PROJECT_NAME})

SET (PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
SET (PROJECT_PYSOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pysource)
SET (PROJECT_LIBS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libs)
SET (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

FILE (GLOB sources ${PROJECT_SOURCE_DIR}/*.cpp)
FILE (GLOB headers ${PROJECT_SOURCE_DIR}/*.h)
FILE (GLOB_RECURSE pythons ${PROJECT_PYSOURCE_DIR}/*.py)

# group python sources
foreach (resource ${pythons})
  get_filename_component(currentdir ${resource} PATH)
  list (APPEND includeDirs ${currentdir})
  get_filename_component(relative ${currentdir} PATH)
  STRING(REGEX REPLACE ${CMAKE_SOURCE_DIR} "" groupname "${currentdir}")
  STRING(REGEX REPLACE "/" "\\\\" groupname "${groupname}")
  source_group(${groupname} FILES ${resource})
endforeach()


# ---- Lars libs ----

ADD_SUBDIRECTORY(${PROJECT_LIBS_DIR}/lars)t.boeddeker
# ---- Finite Differences ----

LIST(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/python.cpp ${PROJECT_SOURCE_DIR}/main.cpp)
add_library(FD STATIC ${sources} ${headers} )

# ---- executable ----

add_executable(run ${PROJECT_SOURCE_DIR}/main.cpp)
TARGET_LINK_LIBRARIES(run FD)

# ---- python ----
SET(Python_ADDITIONAL_VERSIONS 2.7 2.6 2.5)
FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)
FIND_PACKAGE(NumPy)
FIND_PACKAGE(Boost COMPONENTS python)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES} ${PYTHON_LIBRARIES}) # Deprecated but so convenient!

INCLUDE_DIRECTORIES(${PROJECT_LIBS_DIR})


add_definitions(-std=c++11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC"  )
PYTHON_ADD_MODULE(_pypropagate ${PROJECT_SOURCE_DIR}/python.cpp ${pythons})
TARGET_LINK_LIBRARIES(_pypropagate FD)

# ---- install ----
install(TARGETS _pypropagate DESTINATION pypropagate)
# install(FILES ${pythons} DESTINATION pyfinitedifferences_3_2)
install(DIRECTORY ${PROJECT_PYSOURCE_DIR}/ DESTINATION pypropagate)

