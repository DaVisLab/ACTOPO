cmake_minimum_required(VERSION 3.2)

project(ttkCompactTriangulationPreconditioningCmd)

if(TARGET ttkCompactTriangulationPreconditioning)
  add_executable(${PROJECT_NAME} main.cpp)
  target_link_libraries(${PROJECT_NAME}
    PRIVATE
      ttkCompactTriangulationPreconditioning
      VTK::IOXML
    )
  set_target_properties(${PROJECT_NAME}
    PROPERTIES
      INSTALL_RPATH
        "${CMAKE_INSTALL_RPATH}"
    )
  install(
    TARGETS
      ${PROJECT_NAME}
    RUNTIME DESTINATION
      ${TTK_INSTALL_BINARY_DIR}
    )
endif()
