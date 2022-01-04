add_executable(MinimizerTest
	MinimizerTest.cpp
)

target_include_directories(MinimizerTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(MinimizerTest
	Mcs
	${NINJA_MOMRECON_LIBRARIES}
)

install(TARGETS MinimizerTest DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(BetheBlochCheck
	BetheBlochCheck.cpp
)

target_include_directories(BetheBlochCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(BetheBlochCheck
	Mcs
	${NINJA_MOMRECON_LIBRARIES}
)

install(TARGETS BetheBlochCheck DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(DeltaTheta
	DeltaTheta.cpp
)

target_include_directories(DeltaTheta PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DeltaTheta
	Mcs
	${NINJA_MOMRECON_LIBRARIES}
)

install(TARGETS DeltaTheta DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(DeltaThetaLatCheck
	DeltaThetaLatCheck.cpp
)

target_include_directories(DeltaThetaLatCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DeltaThetaLatCheck
	Mcs
	${NINJA_MOMRECON_LIBRARIES}
)

install(TARGETS DeltaThetaLatCheck DESTINATION "${CMAKE_INSTALL_BINDIR}")