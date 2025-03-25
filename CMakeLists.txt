#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(RHICfAnalysis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")

#----------------------------------------------------------------------------
# Find ROOT (required package)
set(ROOTPATH /afs/rhic.bnl.gov/star/ROOT/5.34.38/.sl73_gcc485/rootdeb CACHE PATH "Path to the ROOT library")
message(STATUS "ROOT Path: ${ROOTPATH}")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_PREFIX_PATH ${ROOTPATH})

find_package(ROOT REQUIRED COMPONENTS RIO Net Eve EG RGL Geom Gui Graf Graf3d Ged Matrix MathCore Table Physics Minuit2 Spectrum ROOTDataFrame)

#----------------------------------------------------------------------------
set(SOURCE_DIRECTORY_LIST
  ${CMAKE_CURRENT_LIST_DIR}/source/Base
  ${CMAKE_CURRENT_LIST_DIR}/source/Util

  CACHE INTERNAL ""
)

set(TOOLS_LIBRARIES_LIST ${ROOT_LIBRARIES} CACHE INTERNAL "")

foreach(_directory ${SOURCE_DIRECTORY_LIST})
  include_directories(${_directory})

  string(REPLACE ${CMAKE_SOURCE_DIR} "" tmpString ${_directory})
  string(REPLACE / _ _name ${tmpString})

  list(APPEND SOURCE_NAME_LIST ${_name})
  set(${_name}_G_NAME G__${_name}Dict)
  file(GLOB ${_name}_SOURCE_FILES ${_directory}/*.cc)
  file(GLOB ${_name}_HEADER_FILES ${_directory}/*.hh)
endforeach(_directory)

set(RHICFSIMDST_DIR ${CMAKE_CURRENT_LIST_DIR}/source/StRHICfEventDst)
link_directories(${LINK_DIRECTORIES})

foreach(_directory ${RHICFSIMDST_DIR})
  include_directories(${_directory})

  string(REPLACE ${CMAKE_SOURCE_DIR} "" tmpString ${_directory})
  string(REPLACE / _ _name ${tmpString})

  list(APPEND DIR_NAME_LIST ${_name})
  set(${_name}_G_NAME G__${_name}Dict)
  file(GLOB ${_name}_LINKDEF_FILE ${_directory}/LinkDef.h)
  file(GLOB ${_name}_SOURCE_FILES ${_directory}/St*.cxx)
  file(GLOB ${_name}_HEADER_FILES ${_directory}/St*.h)
endforeach(_directory)

include_directories( ${ROOT_INCLUDE_DIR} )

add_definitions(${ROOT_CXX_FLAGS})

foreach(_name ${DIR_NAME_LIST})
  ROOT_GENERATE_DICTIONARY(${${_name}_G_NAME}
    ${${_name}_HEADER_FILES}
    LINKDEF ${${_name}_LINKDEF_FILE}
  )
endforeach(_name)

foreach(_name ${SOURCE_NAME_LIST})
  foreach(_file ${${_name}_SOURCE_FILES})
    list(APPEND LIB_LIST ${_file})
  endforeach(_file)
endforeach(_name)

foreach(_name ${DIR_NAME_LIST})
  foreach(_file ${${_name}_SOURCE_FILES})
    list(APPEND LIB_LIST ${_file})
  endforeach(_file)
  list(APPEND LIB_LIST ${${_name}_G_NAME}.cxx)
endforeach(_name)


add_library(RHICfAnalysisLib SHARED ${LIB_LIST})
target_link_libraries(RHICfAnalysisLib ${TOOLS_LIBRARIES_LIST})


# =======================  Execute files =========================
add_executable(calculateAN ANCalculation.cc)
target_link_libraries(calculateAN RHICfAnalysisLib)
