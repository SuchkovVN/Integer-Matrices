add_executable(compute "src/compute.cpp")
target_link_libraries(compute PRIVATE intmat::pari intmat::utils)