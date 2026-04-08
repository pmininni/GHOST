add_library(stepper_lib STATIC
  gstepper_base.f90
  canuto_stepper.f90
  gexrk_stepper.f90
  stepper_factory.f90
)

target_include_directories(stepper_lib PUBLIC
  "$<TARGET_PROPERTY:project_fft,INTERFACE_INCLUDE_DIRECTORIES>"
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/mpiio
  ${CMAKE_BINARY_DIR}/src/memmgt
  ${CMAKE_BINARY_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/pdes
  ${CMAKE_BINARY_DIR}/src/particles
  ${CMAKE_BINARY_DIR}/src/${FFTP}
  )

# Make sure dependencies are satisfied
add_dependencies(stepper_lib
  fftp_lib
  pseudo_lib
  mpiio_lib
  memmgt_lib
  pdes_lib
  particles_lib
  )

# Link MPI if needed
target_link_libraries(stepper_lib PUBLIC project_mpi pseudo_lib utils_lib)
