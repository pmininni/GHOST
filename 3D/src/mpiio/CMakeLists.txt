# Try to add binary_io if it exists
#set(_MPIIO_SOURCES)
#if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/binary_io.f90")
#  list(APPEND _MPIIO_SOURCES binary_io.f90)
#endif()

add_library(mpiio_lib STATIC
  binary_mod.f90
  binary_io.f90
)

target_include_directories(mpiio_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/mpiio
  ${CMAKE_BINARY_DIR}/src/${FFTP}
  ${CMAKE_BINARY_DIR}/src/mpiio
  ${CMAKE_BINARY_DIR}/src/utils
  )

# Make sure dependencies are satisfied
add_dependencies(mpiio_lib fftp_lib timer_lib)

# Link MPI if needed
target_link_libraries(mpiio_lib PUBLIC project_mpi)
