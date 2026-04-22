# ---------------------- Analysis tools -------------------------------
# Boots3D:
add_executable(BOOTS "${SRC_DIR}/analysis/boots3D.f90")

add_dependencies(BOOTS fftp_lib pseudo_lib mpiio_lib)

target_link_libraries(BOOTS
    PRIVATE
        project_fft
        project_mpi
        project_includes
        pseudo_lib
        mpiio_lib
)

set_target_properties(BOOTS PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
)
