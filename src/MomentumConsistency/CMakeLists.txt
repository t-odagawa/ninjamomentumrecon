add_executable(MomentumConsistency
	MomentumConsistency.cpp
        )

target_include_directories(MomentumConsistency PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(MomentumConsistency
	Mcs
	Pid
	Range
	${NINJA_MOMRECON_LIBRARIES})

install(TARGETS MomentumConsistency DESTINATION "${CMAKE_INSTALL_BINDIR}")