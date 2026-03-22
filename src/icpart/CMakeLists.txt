add_library(icp_lib STATIC
  icp_base.f90
  icp_position.f90
  #icp_velocity.f90
  icp_factory.f90
)

target_include_directories(icp_lib PUBLIC
  "$<TARGET_PROPERTY:project_fft,INTERFACE_INCLUDE_DIRECTORIES>"
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/mpiio
  ${CMAKE_BINARY_DIR}/src/memmgt
  ${CMAKE_BINARY_DIR}/src/pdes
  ${CMAKE_BINARY_DIR}/src/particles
  ${CMAKE_BINARY_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/utils
  ${CMAKE_BINARY_DIR}/src/${FFTP}
  )

# Make sure dependencies are satisfied
add_dependencies(icp_lib
  fftp_lib
  pseudo_lib
  memmgt_lib
  pdes_lib
  particles_lib
  utils_lib
  )

# Link MPI if needed
target_link_libraries(icp_lib PUBLIC project_mpi utils_lib)
