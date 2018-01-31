cmake_minimum_required (VERSION 3.0)

project (ProteinPipe)
enable_language (CXX)
set (CMAKE_CXX_STANDARD 11)

if(APPLE)
    set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/" ${CMAKE_MODULE_PATH})
endif()

######################################## macro for Visual Studio  ===========================
macro(configure_link_flags)
  set(MSVC_C_CXX_FLAGS
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_MINSIZEREL
    CMAKE_C_FLAGS_RELEASE
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_RELWITHDEBINFO
  )
  if(${VCPKG_TARGET_TRIPLET} MATCHES "static")
    message(STATUS
      "MSVC: static link"
    )
    foreach(flag ${MSVC_C_CXX_FLAGS})
      if(${flag} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag} "${${flag}}")
      endif()
    endforeach()
  else()
    message(STATUS
      "MSVC: dynamic link"
    )
    foreach(flag ${MSVC_C_CXX_FLAGS})
      if(${flag} MATCHES "/MT")
        string(REGEX REPLACE "/MT" "/MD" ${flag} "${${flag}}")
      endif()
    endforeach()
  endif()
endmacro()


if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-DBOOST_NO_RVALUE_REFERENCES)
  add_definitions(-DNOMINMAX)
  configure_link_flags()
endif()

if (APPLE)
	add_definitions(-Wno-deprecated-declarations)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_definitions(-D_POSIX_C_SOURCE)
endif()


# make sure that the default is a RELEASE
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()



### Set output directories on a per-configuration base
# Single configuration
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )
# Multi configuration
foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
    string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_SOURCE_DIR}/bin )
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

###########################################################################################

set (HPP_DIR ${CMAKE_SOURCE_DIR}/cpp/lib) # directory of hpp
set (CPP_DIR ${CMAKE_SOURCE_DIR}/cpp) # cpp directory

find_package(OpenGL REQUIRED)
if(APPLE)
find_package(FREEGLUT REQUIRED)
else()
find_package(GLUT REQUIRED)
endif()
include_directories(${HPP_DIR} ${DATA_ANALYSIS} ${OpenGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR})

if (CYGWIN)
	set (CMAKE_CXX_FLAGS "-std=gnu++14 ${CMAKE_CXX_FLAGS}")
endif()
if (NOT MSVC)
  set (CMAKE_CXX_FLAGS "-ffast-math ${CMAKE_CXX_FLAGS}") # perform more aggressive floating-point optimizations
  set (CMAKE_C_FLAGS "-ffast-math ${CMAKE_C_FLAGS}")
endif()

# ============ .EXE ===============

add_executable(	pdb2xyz
				${CPP_DIR}/pdb2xyz.cpp
			  )
add_executable(	viewer
				${CPP_DIR}/viewer.cpp
				${HPP_DIR}/draw3d.hpp
        ${HPP_DIR}/object.hpp
				${HPP_DIR}/protein.hpp
			  )
target_link_libraries(viewer ${GLUT_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
