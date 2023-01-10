message (STATUS "Test tools...")

add_executable(MomChainTextIO
	MomChainTextIO.cpp
	)


# install the execute in the bin folder
install(TARGETS MomChainTextIO DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(TestPidHist
	TestPidHist.cpp
	)

target_link_libraries(TestPidHist
	Pid
	${NINJA_MOMRECON_LIBRARIES})

# install the execute in the bin folder
install(TARGETS TestPidHist DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(TestPidLine
	TestPidLine.cpp
	)

target_link_libraries(TestPidLine
	Pid
	${NINJA_MOMRECON_LIBRARIES})

# install the execute in the bin folder
install(TARGETS TestPidLine DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(CompareTwoMomch
	CompareTwoMomch.cpp)

target_link_libraries(CompareTwoMomch
	Mcs
	${NINJA_MOMRECON_LIBRARIES})

install(TARGETS CompareTwoMomch DESTINATION ${CMAKE_INSTALL_BINDIR})