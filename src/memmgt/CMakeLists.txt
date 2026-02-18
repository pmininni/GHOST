add_library(memmgt_lib STATIC
  workspace.f90
  gstate_mod.f90
  gpstate_mod.f90)

target_include_directories(memmgt_lib PUBLIC
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/${FFTP}
  )

# Make sure dependencies are satisfied
add_dependencies(memmgt_lib
  fftp_lib
  pseudo_lib
  )

# Link MPI if needed
target_link_libraries(memmgt_lib PUBLIC project_mpi)
