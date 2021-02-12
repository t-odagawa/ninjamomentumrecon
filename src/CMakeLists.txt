add_executable(NinjaMomentumRecon
        main.cpp
        McsCoordinateMethod.hpp
        McsCoordinateMethod.cpp
        McsAngleMethod.hpp
        McsAngleMethod.cpp
        )

target_include_directories(NinjaMomentumRecon PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(NinjaMomentumRecon PRIVATE ${WAGASCI_RECON_LIBRARIES})
install(TARGETS NinjaMomentumRecon DESTINATION "${CMAKE_INSTALL_BINDIR}")