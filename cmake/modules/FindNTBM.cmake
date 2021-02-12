# This module difines
# NTBM_LIBRARY, the name of the library to link against
# NTBM_FOUND, if false, do not try to link to NTBM
# NTBM_INCLUDE_DIR, where to find NTBM headers
#
# Note: If you see and empty NTBM_LIBRARY_TEMP in your ocnfiguration
# and no NTBM_LIBRARY, it means CMake did not find your NTBM library
# (libNTBM,dylib, libNTBM.so etc).
# These values are used to generate the final NTBM_LIBRARY variable,
# but when these values are unset, NTBM_LIBRARY does not get created.
#
#
# ${NINJARECONDIR} is and environment variale that would
# correspond to the CMAKE_INSTALL_PREFIX=${NINJARECONDIR}
# used in buildng NTBM.

SET(NTBM_SEARCH_PATHS
  /usr/local
  /opt/local
  ${NINJARECON_PATH}
  )

FIND_PATH(NTBM_INCLUDE_DIR NTBMSummary.hh
  HINTS
  $ENV{NINJARECONDIR}
  PATH_SUFFIXES include/ninja/recon include
  PATHS ${NTBM_SEARCH_PATHS}
  )

FIND_LIBRARY(NTBM_LIBRARY_TEMP
  NAMES NTBM
  HINTS
  $ENV{NINJARECONDIR}
  PATH_SUFFIXES lib64/ninja/recon lib/ninja/recon
  PATHS ${NTBM_SEARCH_PATHS}
  )

IF(NTBM_LIBRARY_TEMP)
  # Set the final string heere so teh GUI reflects the final state.
  SET(NTBM_LIBRARY ${NTBM_LIBRARY_TEMP} CACHE STRING "Where the NTBM Library can be found")
  # Set the temp bariable to INTERNAL so it is not seen in the CMake GUI
  SET(NTBM_LIBRARY_TEMP "${NTBM_LIBRARY_TEMP}" CACHE INTERNAL "")
ENDIF(NTBM_LIBRARY_TEMP)

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(NTBM REQUIRED_VARS NTBM_LIBRARY NTBM_INCLUDE_DIR)
