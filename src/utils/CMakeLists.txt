add_library(timer_lib STATIC gtimer_mod.f90)
add_library(utils_lib STATIC
  gutils_mod.F90
  gbench_mod.f90
)

# Include directories
target_include_directories(timer_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/utils
  ${CMAKE_BINARY_DIR}
)
target_include_directories(utils_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/utils
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/${FFTP}
)

# Make sure dependencies are satisfied
add_dependencies(utils_lib timer_lib fftp_lib pseudo_lib)

# Link MPI if needed
target_link_libraries(timer_lib PUBLIC project_mpi)
target_link_libraries(utils_lib PUBLIC project_mpi project_fft timer_lib pseudo_lib)
