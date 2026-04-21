add_library(mpiio_lib STATIC
  binary_mod.f90
  binary_io.f90
)

target_include_directories(mpiio_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/mpiio
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/${FFTP}
  ${CMAKE_BINARY_DIR}/src/mpiio
  ${CMAKE_BINARY_DIR}/src/utils
  )

# Make sure dependencies are satisfied
add_dependencies(mpiio_lib fftp_lib utils_lib timer_lib)

# Link MPI if needed
target_link_libraries(mpiio_lib PUBLIC project_mpi)
