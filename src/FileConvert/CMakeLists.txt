add_executable(FileConvert
	FileConvert.cpp
        )

add_executable(FileReConvert
	FileReConvert.cpp
        )

add_executable(FileConvertNonTime
	FileConvertNonTime.cpp
	)

add_executable(FileReConvertNonTime
	FileReConvertNonTime.cpp
	)

add_executable(DailySeparateMomChain
	DailySeparateMomChain.cpp
        )

add_executable(DailyMergeMomChain
	DailyMergeMomChain.cpp
	)

add_executable(MCOutputToChain
	MCOutputToChain.cpp
	)

add_executable(AddBabyMindToChain
	AddBabyMindToChain.cpp)

add_executable(NuFactEmergency
	NuFactEmergency.cpp)

add_executable(NuFactHistogram
	NuFactHistogram.cpp)

add_executable(AddBabyMindToChainMcBase
	AddBabyMindToChainMcBase.cpp)
	
add_executable(PhdThesisHistogram
	PhdThesisHistogram.cpp)

add_executable(PhdThesisHistogram0pi
	PhdThesisHistogram0pi.cpp)

add_executable(Phd0pi1p
	Phd0pi1p.cpp)

add_executable(Phd0pi2p
	Phd0pi2p.cpp)

add_executable(Phd1pi0p
	Phd1pi0p.cpp)

target_include_directories(FileConvert PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(FileReConvert PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(FileConvertNonTime PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(FileReConvertNonTime PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(DailySeparateMomChain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(DailyMergeMomChain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(MCOutputToChain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(AddBabyMindToChain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(NuFactEmergency PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})					    
target_include_directories(NuFactHistogram PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(AddBabyMindToChainMcBase PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(PhdThesisHistogram PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(PhdThesisHistogram0pi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Phd0pi1p PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Phd0pi2p PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Phd1pi0p PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(FileConvert
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(FileReConvert
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(FileConvertNonTime
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(FileReConvertNonTime
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(DailySeparateMomChain
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(DailyMergeMomChain
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(MCOutputToChain
	Mcs
	Connection
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(AddBabyMindToChain
	Mcs
	Match
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(NuFactEmergency
	Mcs
	Range
	Match
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(NuFactHistogram
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(AddBabyMindToChainMcBase
	Mcs
	Match
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(PhdThesisHistogram
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(PhdThesisHistogram0pi
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(Phd0pi1p
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(Phd0pi2p
	Mcs
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(Phd1pi0p
	Mcs
	${NINJA_MOMRECON_LIBRARIES})

install(TARGETS FileConvert DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS FileReConvert DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS FileConvertNonTime DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS FileReConvertNonTime DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS DailySeparateMomChain DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS DailyMergeMomChain DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS MCOutputToChain DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS AddBabyMindToChain DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS NuFactEmergency DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS NuFactHistogram DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS AddBabyMindToChainMcBase DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS PhdThesisHistogram DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS PhdThesisHistogram0pi DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS Phd0pi1p DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS Phd0pi2p DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS Phd1pi0p DESTINATION "${CMAKE_INSTALL_BINDIR}")
