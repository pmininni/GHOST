add_library(pseudo_lib STATIC
  pseudospec3D_mod.f90
  pseudospec3D_hd.f90
  pseudospec3D_norm.f90
  pseudospec3D_phd.f90
  pseudospec3D_rot.f90
  pseudospec3D_mhd.f90
)

target_include_directories(pseudo_lib PUBLIC
  "$<TARGET_PROPERTY:project_fft,INTERFACE_INCLUDE_DIRECTORIES>"
   ${CMAKE_BINARY_DIR}
   ${CMAKE_BINARY_DIR}/src/utils
   ${CMAKE_BINARY_DIR}/src/${FFTP}
  )

# Make sure dependencies are satisfied
add_dependencies(pseudo_lib
  fftp_lib
  )

# Link MPI if needed
target_link_libraries(pseudo_lib PUBLIC project_mpi fftp_lib utils_lib)
