add_executable(FileConvert
	FileConvert.cpp
#	FileConvert.hpp
        )

target_include_directories(FileConvert PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(FileConvert
	${ROOT_LIBRARIES}
	${NINJA_MOMRECON_LIBRARIES})
install(TARGETS FileConvert DESTINATION "${CMAKE_INSTALL_BINDIR}")