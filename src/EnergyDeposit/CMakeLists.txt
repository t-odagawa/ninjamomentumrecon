add_executable(EnergyDeposit
	EnergyDeposit.cpp
#	EnergyDeposit.hpp
        )

target_include_directories(EnergyDeposit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(EnergyDeposit
	${ROOT_LIBRARIES}
	${NINJA_MOMRECON_LIBRARIES})
install(TARGETS EnergyDeposit DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(MinimizerTest
	MinimizerTest.cpp
	MinimizerTest.hpp	
)

target_include_directories(MinimizerTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(MinimizerTest
#	${ROOT_LIBRARIES}
	Minuit
	${NINJA_MOMRECON_LIBRARIES})

install(TARGETS MinimizerTest DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(BetheBlochCheck
	BetheBlochCheck.cpp
)

target_include_directories(BetheBlochCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(BetheBlochCheck
	${NINJA_MOMRECON_LIBRARIES}
)

install(TARGETS BetheBlochCheck DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(DeltaTheta
	DeltaTheta.cpp
)

target_include_directories(DeltaTheta PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DeltaTheta
	${NINJA_MOMRECON_LIBRARIES}
)

install(TARGETS DeltaTheta DESTINATION "${CMAKE_INSTALL_BINDIR}")