# This module difines
# MCS_LIBRARY, the name of the library to link against
# MCS_FOUND, if false, do not try to link to NTBM
# MCS_INCLUDE_DIR, where to find NTBM headers
#
# Note: If you see and empty MCS_LIBRARY_TEMP in your ocnfiguration
# and no MCS_LIBRARY, it means CMake did not find your MCS library
# (libMCS,dylib, libMCS.so etc).
# These values are used to generate the final MCS_LIBRARY variable,
# but when these values are unset, MCS_LIBRARY does not get created.
#
#
# ${NINJARECONDIR} is and environment variale that would
# correspond to the CMAKE_INSTALL_PREFIX=${NINJARECONDIR}
# used in buildng MCS.

SET(MCS_SEARCH_PATHS
  /usr/local
  /opt/local
  ${NINJARECON_PATH}
  )

FIND_PATH(MCS_INCLUDE_DIR MCSSummary.hh
  HINTS
  $ENV{NINJARECONDIR}
  PATH_SUFFIXES include/ninja/recon include
  PATHS ${NTBM_SEARCH_PATHS}
  )

FIND_LIBRARY(MCS_LIBRARY_TEMP
  NAMES MCS
  HINTS
  $ENV{NINJARECONDIR}
  PATH_SUFFIXES lib64/ninja/recon lib/ninja/recon
  PATHS ${NTBM_SEARCH_PATHS}
  )

IF(MCS_LIBRARY_TEMP)
  # Set the final string heere so teh GUI reflects the final state.
  SET(MCS_LIBRARY ${MCS_LIBRARY_TEMP} CACHE STRING "Where the MCS Library can be found")
  # Set the temp bariable to INTERNAL so it is not seen in the CMake GUI
  SET(MCS_LIBRARY_TEMP "${MCS_LIBRARY_TEMP}" CACHE INTERNAL "")
ENDIF(MCS_LIBRARY_TEMP)

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MCS REQUIRED_VARS MCS_LIBRARY MCS_INCLUDE_DIR)
