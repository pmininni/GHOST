# ---------------------- GHOST tools ----------------------------------
# Boots3D:
add_executable(BOOTS "${SRC_DIR}/tools/boots3D.f90")

# pseudo_lib comes in through utils_lib and mpiio_lib; listing it here
# too puts it twice on the link line (ld warns about the duplicate).
target_link_libraries(BOOTS
    PRIVATE
        project_includes
        project_gpu
        utils_lib
        mpiio_lib
)

set_target_properties(BOOTS PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
)

# Python interface:
add_library(ghost_python SHARED "${SRC_DIR}/tools/gpython_lib.f90")

set_target_properties(ghost_python PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${LIB_DIR}"
)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_dependencies(ghost_python stepper_lib)

target_link_libraries(ghost_python
    PRIVATE
        project_includes
        project_gpu
        mpiio_lib
        particles_lib
        ic_lib
        force_lib
        icp_lib
        stepper_lib
)
