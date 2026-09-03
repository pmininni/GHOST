# pseudobase: the basic modules (grid, wavenumbers, box, fft plans,
# random numbers), needed by the memory management before the
# pseudospectral kernels can be built
add_library(pseudobase_lib STATIC
  pseudospec3D_mod.f90
)

target_include_directories(pseudobase_lib PUBLIC
  "$<TARGET_PROPERTY:project_fft,INTERFACE_INCLUDE_DIRECTORIES>"
   ${CMAKE_BINARY_DIR}
   ${CMAKE_BINARY_DIR}/src/utils
   ${CMAKE_BINARY_DIR}/src/${FFTP}
   ${CMAKE_BINARY_DIR}/src/memmgt
  )

add_dependencies(pseudobase_lib fftp_lib gmem_lib)
target_link_libraries(pseudobase_lib PUBLIC project_mpi fftp_lib gmem_lib)

# pseudo: the pseudospectral kernels; they take their temporaries
# from the workspace pool (memmgt_lib)
add_library(pseudo_lib STATIC
  pseudospec3D_hd.f90
  pseudospec3D_norm.f90
  pseudospec3D_phd.f90
  pseudospec3D_rot.f90
  pseudospec3D_prot.f90
  pseudospec3D_mhd.f90
  pseudospec3D_bouss.f90
)

target_include_directories(pseudo_lib PUBLIC
  "$<TARGET_PROPERTY:project_fft,INTERFACE_INCLUDE_DIRECTORIES>"
   ${CMAKE_BINARY_DIR}
   ${CMAKE_BINARY_DIR}/src/utils
   ${CMAKE_BINARY_DIR}/src/${FFTP}
   ${CMAKE_BINARY_DIR}/src/memmgt
  )

# Make sure dependencies are satisfied
add_dependencies(pseudo_lib
  fftp_lib
  pseudobase_lib
  memmgt_lib
  )

# Link MPI if needed
target_link_libraries(pseudo_lib PUBLIC project_mpi fftp_lib pseudobase_lib memmgt_lib)
