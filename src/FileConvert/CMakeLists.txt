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