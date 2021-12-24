message (STATUS "Test tools...")

add_executable(MomChainTextIO
	MomChainTextIO.cpp
	)


# install the execute in the bin folder
install(TARGETS MomChainTextIO DESTINATION ${CMAKE_INSTALL_BINDIR})