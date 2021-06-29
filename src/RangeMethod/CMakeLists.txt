add_executable(RangeEnergyDepositCheck
	RangeEnergyDepositCheck.cpp
        )

target_include_directories(RangeEnergyDepositCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(RangeEnergyDepositCheck
	${ROOT_LIBRARIES}
	${NINJA_MOMRECON_LIBRARIES})
install(TARGETS RangeEnergyDepositCheck DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(TestReadOut
	TestReadOut.cpp
        )

target_include_directories(TestReadOut PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TestReadOut
	${ROOT_LIBRARIES}
	${NINJA_MOMRECON_LIBRARIES})
install(TARGETS TestReadOut DESTINATION "${CMAKE_INSTALL_BINDIR}")
