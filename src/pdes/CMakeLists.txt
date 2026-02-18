add_library(pdes_lib STATIC
  equation_base.f90
  hd_solver.f90
  mhd_solver.f90
  equation_factory.f90
)

target_include_directories(pdes_lib PUBLIC
  "$<TARGET_PROPERTY:project_fft,INTERFACE_INCLUDE_DIRECTORIES>"
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/mpiio
  ${CMAKE_BINARY_DIR}/src/memmgt
  ${CMAKE_BINARY_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/${FFTP}
  )

# Make sure dependencies are satisfied
add_dependencies(pdes_lib
  fftp_lib
  pseudo_lib
  mpiio_lib
  memmgt_lib
  )

# Link MPI if needed
target_link_libraries(pdes_lib PUBLIC project_mpi pseudo_lib utils_lib)
