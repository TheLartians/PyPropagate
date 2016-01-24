
file (GLOB headers "*.h")

add_library(LarsHeaderLibs EXCLUDE_FROM_ALL ${headers})
SET_TARGET_PROPERTIES(LarsHeaderLibs PROPERTIES LINKER_LANGUAGE CXX)

