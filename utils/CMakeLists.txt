set(TARGET_NAME "intmat_utils")

file(GLOB_RECURSE SRCS "*.h" "*.cpp")
add_library(${TARGET_NAME} STATIC ${SRCS})
add_library(intmat::utils ALIAS ${TARGET_NAME})

target_include_directories(${TARGET_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
