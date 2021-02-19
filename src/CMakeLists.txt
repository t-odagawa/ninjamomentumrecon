add_executable(NinjaMomentumRecon
        McsCoordinateMethod.hpp
        McsCoordinateMethod.cpp
        McsAngleMethod.hpp
        McsAngleMethod.cpp
	McsCommon.hpp
	McsCommon.cpp
        main.cpp
        )

target_include_directories(NinjaMomentumRecon PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(NinjaMomentumRecon
	${ROOT_LIBRARIES}
	${NINJA_MOMRECON_LIBRARIES})
install(TARGETS NinjaMomentumRecon DESTINATION "${CMAKE_INSTALL_BINDIR}")