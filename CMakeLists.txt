
#project
#version salome version 9.3

# this cmake file comes from FreeCAD
find_package(OpenCasCade REQUIRED)
include_directories(
  ${OCC_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

# if other occt detection cmake approaches is used, please define  ${OCC_LIBRARY_DIR} ${OCC_LIBRARIES}
# for subdirectory cmakefiles.txt

add_subdirectory("GEOMAlgo")