# Path to generated Fortran files
set(GEN_PRECIS ${CMAKE_BINARY_DIR}/fprecision_mod.f90)
set(GEN_FFTMOD ${CMAKE_BINARY_DIR}/fftp_mod.f90)
set(GEN_FFTP3D ${CMAKE_BINARY_DIR}/fftp3D.f90)

add_library(fftp_lib STATIC
  fftw_wrappers.h
  fftp.f90
  ${GEN_PRECIS}
  ${GEN_FFTMOD}
  ${GEN_FFTP3D}
)

# Mark generated files
set_source_files_properties(${GEN_PRECIS} PROPERTIES GENERATED TRUE)
set_source_files_properties(${GEN_FFTMOD} PROPERTIES GENERATED TRUE)
set_source_files_properties(${GEN_FFTP3D} PROPERTIES GENERATED TRUE)

# Include directories
target_include_directories(fftp_lib PUBLIC
  "$<TARGET_PROPERTY:project_fft,INTERFACE_INCLUDE_DIRECTORIES>"
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src/utils
)

# Make sure dependencies are satisfied
add_dependencies(fftp_lib project_preprocess timer_lib)

# Link MPI if needed
target_link_libraries(fftp_lib PUBLIC project_mpi)
