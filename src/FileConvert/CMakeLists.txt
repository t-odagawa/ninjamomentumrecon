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

target_include_directories(FileConvert PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(FileReConvert PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(FileConvertNonTime PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(FileReConvertNonTime PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(DailySeparateMomChain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(DailyMergeMomChain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

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
install(TARGETS FileConvert DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS FileReConvert DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS FileConvertNonTime DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS FileReConvertNonTime DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS DailySeparateMomChain DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS DailyMergeMomChain DESTINATION "${CMAKE_INSTALL_BINDIR}")