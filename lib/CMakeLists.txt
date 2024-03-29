add_library(Mcs SHARED
  McsConst.hpp
  McsClass.hpp
  McsClass.cpp
  McsFunction.hpp
  McsFunction.cpp
)

target_include_directories(Mcs PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB MCS_LIB_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
set_target_properties(Mcs PROPERTIES PUBLIC_HEADER "${MCS_LIB_INCLUDES}")

install(TARGETS Mcs
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCDIR}")

add_library(Range SHARED
		  RangeSpline.hpp
		  RangeSpline.cpp
		  RangeFunction.hpp
		  RangeFunction.cpp
)

target_include_directories(Range PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB RANGE_LIB_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
set_target_properties(Range PROPERTIES PUBLIC_HEADER "${RANGE_LIB_INCLUDES}")

add_library(Connection SHARED
		       ConnectionData.hpp
		       ConnectionData.cpp
		       ConnectionClass.hpp
		       ConnectionClass.cpp
		       ConnectionFunction.hpp
		       ConnectionFunction.cpp
)

target_include_directories(Connection PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB CONNECTION_LIB_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
set_target_properties(Connection PROPERTIES PUBLIC_HEADER "${CONNECTION_LIB_INCLUDES}")

install(TARGETS Connection
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCDIR}")

add_library(Pid SHARED
		PidData.hpp
		PidData.cpp
		PidClass.hpp
		PidClass.cpp
		PidFunction.hpp
		PidFunction.cpp
)

target_include_directories(Pid PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB PID_LIB_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
set_target_properties(Pid PROPERTIES PUBLIC_HEADER "${PID_LIB_INCLUDES}")

install(TARGETS Pid
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCDIR}")

add_library(Match SHARED
		  MatchData.hpp
		  MatchData.cpp
		  MatchFunction.hpp
		  MatchFunction.cpp
)

target_include_directories(Match PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB MATCH_LIB_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
set_target_properties(Match PROPERTIES PUBLIC_HEADER "${MATCH_LIB_INCLUDES}")

install(TARGETS Match
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCDIR}")