# This module defines
# B2MC_LIBRARY, the name of the library to link against
# B2MC_FOUND, if false, do not try to link to B2MC
# B2MC_INCLUDE_DIR, where to find B2MC headers
#
# Note: If you see an empty B2MC_LIBRARY_TEMP in your configuration
# and no B2MC_LIBRARY, it means CMake did not find your B2MC library
# (libB2MC.dylib, libB2MC.so, etc).
# These values are used to generate the final B2MC_LIBRARY variable,
# but when these values are unset, B2MC_LIBRARY does not get created.
#
#
# $B2MCDIR is an environment variable that would
# correspond to the CMAKE_INSTALL_PREFIX=$B2MCDIR
# used in building B2MC.

SET(B2MC_SEARCH_PATHS
        /usr/local
        /opt/local
        ${B2MC_PATH}
        )

FIND_PATH(B2MC_INCLUDE_DIR B2SpillSummary.hh
        HINTS
        $ENV{B2MCDIR}
        PATH_SUFFIXES include/wagasci/b2 include
        PATHS ${B2MC_SEARCH_PATHS}
        )

FIND_LIBRARY(B2MC_LIBRARY_TEMP
        NAMES B2MC
        HINTS
        $ENV{B2MCDIR}
        PATH_SUFFIXES lib64/wagasci/b2 lib/wagasci/b2
        PATHS ${B2MC_SEARCH_PATHS}
        )

IF(B2MC_LIBRARY_TEMP)
    # Set the final string here so the GUI reflects the final state.
    SET(B2MC_LIBRARY ${B2MC_LIBRARY_TEMP} CACHE STRING "Where the B2MC Library can be found")
    # Set the temp variable to INTERNAL so it is not seen in the CMake GUI
    SET(B2MC_LIBRARY_TEMP "${B2MC_LIBRARY_TEMP}" CACHE INTERNAL "")
ENDIF(B2MC_LIBRARY_TEMP)

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(B2MC REQUIRED_VARS B2MC_LIBRARY B2MC_INCLUDE_DIR)