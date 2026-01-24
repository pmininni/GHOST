add_library(force_lib STATIC
  force_base.f90
  force_velocity.f90
  force_factory.f90
)

target_include_directories(force_lib PUBLIC
  "$<TARGET_PROPERTY:project_fft,INTERFACE_INCLUDE_DIRECTORIES>"
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/mpiio
  ${CMAKE_BINARY_DIR}/src/memmgt
  ${CMAKE_BINARY_DIR}/src/pdes
  ${CMAKE_BINARY_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/utils
  ${CMAKE_BINARY_DIR}/src/${FFTP}
  )

# Make sure dependencies are satisfied
add_dependencies(force_lib
  fftp_lib
  pseudo_lib
  memmgt_lib
  pdes_lib
  utils_lib
  )

# Link MPI if needed
target_link_libraries(force_lib PUBLIC project_mpi utils_lib)
