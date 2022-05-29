add_executable(PidRange
	PidRange.cpp
        )

target_include_directories(PidRange PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(PidRange	
	Mcs
	Pid
	Range
	${NINJA_MOMRECON_LIBRARIES})

install(TARGETS PidRange DESTINATION "${CMAKE_INSTALL_BINDIR}")