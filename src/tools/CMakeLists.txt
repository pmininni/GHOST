# ---------------------- GHOST tools ----------------------------------
# Boots3D:
add_executable(BOOTS "${SRC_DIR}/tools/boots3D.f90")

# pseudo_lib comes in through utils_lib and mpiio_lib
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

# The order matters for the static archives: stepper_lib references the
# ParticleBase and EquationBase classes, so particles_lib and the other
# archives that provide them must follow it (otherwise the shared
# library is left with an undefined virtual table and fails to load).
# Each archive must appear once: with amdflang a duplicate makes the
# device link fail with duplicate symbols.
target_link_libraries(ghost_python
    PRIVATE
        project_includes
        project_gpu
        stepper_lib
        icp_lib
        ic_lib
        force_lib
        particles_lib
        mpiio_lib
)
