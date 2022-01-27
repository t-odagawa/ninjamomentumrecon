add_executable(RangeEnergyDepositCheck
	RangeEnergyDepositCheck.cpp
        )

target_include_directories(RangeEnergyDepositCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(RangeEnergyDepositCheck
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
install(TARGETS RangeEnergyDepositCheck DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(RangeMethod
	RangeMethod.cpp
	)

target_include_directories(RangeMethod PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(RangeMethod
	Mcs
	Range
	${NINJA_MOMRECON_LIBRARIES})
install(TARGETS RangeMethod DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(TestReadOut
	TestReadOut.cpp
        )

target_include_directories(TestReadOut PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TestReadOut
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
install(TARGETS TestReadOut DESTINATION "${CMAKE_INSTALL_BINDIR}")
