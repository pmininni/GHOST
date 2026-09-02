add_library(particles_lib STATIC
  pstatus_mod.f90
  gpdbrec_mod.f90
  gpartcomm_mod.f90
  gpsplineint_mod.F90
  particle_base.f90
  lagpart_solver.f90
  inerpart_solver.f90
  particle_factory.f90
)

target_include_directories(particles_lib PUBLIC
  "$<TARGET_PROPERTY:project_fft,INTERFACE_INCLUDE_DIRECTORIES>"
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/memmgt
  ${CMAKE_BINARY_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/pdes
  ${CMAKE_BINARY_DIR}/src/${FFTP}
  )

# Make sure dependencies are satisfied
add_dependencies(particles_lib
  pdes_lib
  )

# Link MPI if needed
target_link_libraries(particles_lib PUBLIC project_mpi pseudo_lib utils_lib pdes_lib)
