add_library(pdes_lib STATIC
  gstate_mod.f90
  equation_base.f90
  hd_solver.f90
)

target_include_directories(pdes_lib PUBLIC
  "$<TARGET_PROPERTY:project_fft,INTERFACE_INCLUDE_DIRECTORIES>"
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/memmgt
  ${CMAKE_BINARY_DIR}/src/${FFTP}
  )

# Make sure dependencies are satisfied
add_dependencies(pdes_lib
  fftp_lib
  memmgt_lib
  pseudo_lib
  )

# Link MPI if needed
target_link_libraries(pdes_lib PUBLIC project_mpi)
