add_library(ic_lib STATIC
  ic_base.f90
  ic_velocity.f90
  ic_factory.f90
)

target_include_directories(ic_lib PUBLIC
  "$<TARGET_PROPERTY:project_fft,INTERFACE_INCLUDE_DIRECTORIES>"
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/mpiio
  ${CMAKE_BINARY_DIR}/src/memmgt
  ${CMAKE_BINARY_DIR}/src/pdes
  ${CMAKE_BINARY_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/${FFTP}
  )

# Make sure dependencies are satisfied
add_dependencies(ic_lib
  fftp_lib
  pseudo_lib
  memmgt_lib
  pdes_lib
  )

# Link MPI if needed
target_link_libraries(ic_lib PUBLIC project_mpi)
