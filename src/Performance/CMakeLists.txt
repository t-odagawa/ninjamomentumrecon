add_executable(BMRangeEccMcs
	BMRangeEccMcs.cpp
        )

target_include_directories(BMRangeEccMcs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(BMRangeEccMcs
	Mcs
	${NINJA_MOMRECON_LIBRARIES})

install(TARGETS BMRangeEccMcs DESTINATION "${CMAKE_INSTALL_BINDIR}")