add_executable(NinjaMomentumRecon
        McsCoordinateMethod.hpp
        McsCoordinateMethod.cpp
        McsAngleMethod.hpp
        McsAngleMethod.cpp
        main.cpp
        )

target_include_directories(NinjaMomentumRecon PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(NinjaMomentumRecon
	${NINJA_MOMRECON_LIBRARIES})
install(TARGETS NinjaMomentumRecon DESTINATION "${CMAKE_INSTALL_BINDIR}")