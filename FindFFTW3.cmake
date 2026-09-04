# FindFFTW3.cmake
# Helper module to locate FFTW3 with optional precision.
#
# Provided variables on success:
#   FFTW3_FOUND
#   FFTW3_INCLUDE_DIRS
#   FFTW3_LIBRARIES        # double precision (fftw3)
#   FFTW3F_LIBRARIES       # single precision (fftw3f)
#   FFTW3_LIBRARY          # raw double-precision library path
#   FFTW3F_LIBRARY         # raw single-precision library path
#
# Hints:
#   Set FFTW3_ROOT to help the search.
#   Example:
#     -DFFTW3_ROOT=/opt/fftw/3.3.10

include_guard(GLOBAL)

# Allow multiple common hint variables
set(_FFTW3_HINTS "")
foreach(hint_var FFTW3_ROOT)
  if(DEFINED ${hint_var})
    list(APPEND _FFTW3_HINTS "${${hint_var}}")
  endif()
endforeach()

# Search for header
find_path(FFTW3_INCLUDE_DIR
  NAMES fftw3.h
  HINTS ${_FFTW3_HINTS}
  PATH_SUFFIXES include
)

# Search for double-precision library
find_library(FFTW3_LIBRARY
  NAMES fftw3
  HINTS ${_FFTW3_HINTS}
  PATH_SUFFIXES lib lib64
)

# Search for single-precision library
find_library(FFTW3F_LIBRARY
  NAMES fftw3f
  HINTS ${_FFTW3_HINTS}
  PATH_SUFFIXES lib lib64
)

# Compose convenience variables
set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
set(FFTW3F_LIBRARIES ${FFTW3F_LIBRARY})

include(FindPackageHandleStandardArgs)
# We consider the package found if we have headers and at least one of the libraries.
find_package_handle_standard_args(FFTW3
  REQUIRED_VARS FFTW3_INCLUDE_DIR
  # Accept either double or single precision library presence
  HANDLE_COMPONENTS
)

# If headers found but neither library was found, mark NOTFOUND
if(FFTW3_INCLUDE_DIR AND NOT FFTW3_LIBRARY AND NOT FFTW3F_LIBRARY)
  set(FFTW3_FOUND FALSE)
endif()

# Mark as advanced
mark_as_advanced(
  FFTW3_INCLUDE_DIR
  FFTW3_LIBRARY
  FFTW3F_LIBRARY
)

# Usage notes:
# - For double precision, link with ${FFTW3_LIBRARIES}
# - For single precision, link with ${FFTW3F_LIBRARIES}
# - Add include dirs: target_include_directories(tgt PUBLIC ${FFTW3_INCLUDE_DIRS})
