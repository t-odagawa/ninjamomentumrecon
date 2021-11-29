add_executable(FileConvert
	FileConvert.cpp
	FileConvert.hpp
        )

add_executable(FileConvertNonTime
	FileConvertNonTime.cpp
	)

add_executable(DailySeparateMomChain
	DailySeparateMomChain.cpp
        )

target_include_directories(FileConvert PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(FileConvertNonTime PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(DailySeparateMomChain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(FileConvert
	${ROOT_LIBRARIES}
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(FileConvertNonTime
	${ROOT_LIBRARIES}
	${NINJA_MOMRECON_LIBRARIES})
target_link_libraries(DailySeparateMomChain
	${ROOT_LIBRARIES}
	${NINJA_MOMRECON_LIBRARIES})

install(TARGETS FileConvert DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS FileConvertNonTime DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS DailySeparateMomChain DESTINATION "${CMAKE_INSTALL_BINDIR}")