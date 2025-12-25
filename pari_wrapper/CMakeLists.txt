set(TARGET_NAME pari_wrapper)

file(GLOB_RECURSE SRCS "*.h" "*.cpp")

add_library(${TARGET_NAME} STATIC ${SRCS})
add_library(intmat::pari ALIAS ${TARGET_NAME})
target_include_directories(${TARGET_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)


target_link_libraries(${TARGET_NAME} PUBLIC ${PARI_LIB})

add_executable(pari_test_exe src/pari_test.cpp)
target_link_libraries(pari_test_exe PUBLIC ${TARGET_NAME})
