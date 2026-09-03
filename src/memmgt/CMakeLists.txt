# gdevice: device management, no dependencies inside GHOST
add_library(gdevice_lib STATIC gdevice_mod.f90)
target_include_directories(gdevice_lib PUBLIC ${CMAKE_BINARY_DIR})
target_link_libraries(gdevice_lib PUBLIC project_gpu)

# gmem: typed allocator, needs the precision kind from fftp_lib
add_library(gmem_lib STATIC gmem_mod.f90)
target_include_directories(gmem_lib PUBLIC
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/${FFTP}
  )
add_dependencies(gmem_lib gdevice_lib fftp_lib)
target_link_libraries(gmem_lib PUBLIC gdevice_lib fftp_lib)

# memmgt: states and workspace pool
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
  gmem_lib
  pseudobase_lib
  )

# Link MPI if needed
target_link_libraries(memmgt_lib PUBLIC project_mpi gmem_lib pseudobase_lib)
