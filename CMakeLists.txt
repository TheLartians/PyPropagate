
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

# ---- Lars libs ----

ADD_SUBDIRECTORY(${PROJECT_LIBS_DIR}/lars)

# ---- python ----
SET(Python_ADDITIONAL_VERSIONS 2.7 2.6 2.5)
FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)
FIND_PACKAGE(NumPy)
FIND_PACKAGE(Boost COMPONENTS python)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES} ${PYTHON_LIBRARIES}) # Deprecated but so convenient!

INCLUDE_DIRECTORIES(${PROJECT_LIBS_DIR})

LIST(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/main.cpp)
LIST(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/python.cpp)

add_definitions(-std=c++11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC"  )
add_library(FD STATIC ${sources} ${headers} )

PYTHON_ADD_MODULE(_finitedifferences ${PROJECT_SOURCE_DIR}/python.cpp ${pythons})
add_executable(FiniteDifferences ${PROJECT_SOURCE_DIR}/main.cpp)

TARGET_LINK_LIBRARIES(_finitedifferences FD)
TARGET_LINK_LIBRARIES(FiniteDifferences FD)

# ---- install ----
install(TARGETS _finitedifferences FD DESTINATION pyfinitedifferences_3_2)
# install(FILES ${pythons} DESTINATION pyfinitedifferences_3_2)
install(DIRECTORY ${PROJECT_PYSOURCE_DIR}/ DESTINATION pyfinitedifferences_3_2)

