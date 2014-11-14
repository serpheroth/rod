project(rod)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.9)

set(ENABLE_CUDA "false")
#set(ENABLE_CUDA "true")
# SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -I${CMAKE_SOURCE_DIR}/lib)
IF(WIN32)
  include_directories("c:\\lib\\include")
  link_directories("c:\\lib\\lib\\x64")
  SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "C:\\Qt\\Qt5.2.1\\5.2.1\\msvc2012_64_opengl")
  SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "C:\\lib")
  add_definitions(-DCURRENT_DIRECTORY=${CMAKE_CURRENT_SOURCE_DIR}/data/)
  set(CMAKE_MAKE_PROGRAM "${CMAKE_MAKE_PROGRAM} -j 8")
# add_definitions(-DDATA_DIRECTORY=\"${CMAKE_CURRENT_SOURCE_DIR}/\" )
  set(CORE_LIBS gslcblas GlU32 freeglut glew32 OpenGL32)
ELSE()
  set(CORE_LIBS  gslcblas)
  add_definitions(-DCURRENT_DIRECTORY=\"${CMAKE_CURRENT_SOURCE_DIR}/\" )
  add_definitions(-DDATA_DIRECTORY=\"${CMAKE_CURRENT_SOURCE_DIR}/data/\" )
  FIND_PACKAGE(OpenGL)
  FIND_PACKAGE(GLEW)
  FIND_PACKAGE(GLUT REQUIRED)
  IF(APPLE)
    SET(CMAKE_CXX_COMPILER "/opt/local/bin/g++-mp-4.7")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -ccbin clang")
    LINK_DIRECTORIES(/opt/local/lib)
    include_directories(/opt/local/include)
    SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "/Users/dj/Qt5.3.2/5.3/clang_64")
  ELSEIF(UNIX)
  ENDIF()
ENDIF()

#------------------------------------------------------------------------------
# Qt
#------------------------------------------------------------------------------
FIND_PACKAGE(Qt5Widgets)
FIND_PACKAGE(Qt5OpenGL)
SET(CMAKE_AUTOMOC ON)
QT5_WRAP_UI(UIS_HDRS src/ui/main_window.ui)
#INCLUDE_DIRECTORIES(${Qt5Widgets_INCLUDE_DIRS})
#INCLUDE_DIRECTORIES(${Qt5OpenGL_INCLUDE_DIRS})
#ADD_DEFINITIONS(${Qt5Widgets_DEFINITIONS})
#ADD_DEFINITIONS(-DQT_PLUGIN)

#---------------------------
# core project
#---------------------------
INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${GLEW_INCLUDE_DIRS})
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include_directories(src/lib/)
include_directories(src/)
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
include_directories(${CMAKE_CURRENT_BINARY_DIR})


#------------------------------------------------------------------------------
# cuda
#------------------------------------------------------------------------------
if(${ENABLE_CUDA} STREQUAL "true")
  FIND_PACKAGE(CUDA REQUIRED)
  message(STATUS "CUDA ENABLED")
  if(WIN32)
    set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "c:\\cuda")
  elseif(APPLE)
    set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "/Developer/NVIDIA/CUDA-6.5")
  #  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-ccbin /opt/local/bin/g++-mp-4.7")
  elseif(UNIX)
  endif()
  add_definitions(-DCUDA_ENABLED)
  #message(STATUS "CUDA_LIBS = " ${CUDA_LIBRARIES})
  #message(STATUS "CUDA_INCLUDES = " ${CUDA_INCLUDE_DIRS})
  #message(STATUS "CUDA SDK = " ${CUDA_TOOLKIT_ROOT_DIR})
  #message(STATUS "CUDA_INCLUDES = " ${CUDA_INCLUDE_DIRS})
  #set(CORE_LIBS  glut GLU GL pthread cudart GLEW)
  include_directories(${CUDA_ROOT_DIRECTORY}/include)
  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_20,code=sm_20")
  SET(CUDA_INCLUDE_DIRS ${CUDA_INCLUDE_DIRS} -I${CMAKE_SOURCE_DIR}/src/lib -I${CMAKE_SOURCE_DIR}/src)
  file(GLOB_RECURSE CUDA_FILES src/*.cu)
  # CUDA_INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
  # message(status  "${CMAKE_SOURCE_DIR}/cuda)
  link_directories(${CUDA_ROOT_DIRECTORY}/lib64)
#  set(CUDA_NVCC_FLAGS_DEBUG "${CUDA_NVCC_FLAGS_DEBUG} -G")
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -G")
  #message(STATUS ${CUDA_NVCC_FLAGS})
  cuda_compile(CUDA_OBJ_FILES ${CUDA_FILES})
else()
  set(CUDA_OBJ_FILES "")
endif()

#------------------------------------------------------------------------------
# cpp files
#------------------------------------------------------------------------------
if(UNIX)
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)
endif()

file(GLOB_RECURSE HEADER_FILES
  src/*.h
)
file(GLOB_RECURSE cpp_files
  src/*.cpp
)
#message(STATUS "headers=" ${HEADER_FILES})
#message(STATUS "cpp files=" ${cpp_files})

find_package(OpenMP)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

ADD_EXECUTABLE(
  ${PROJECT_NAME} 
  ${HEADER_FILES}
  ${UIS_HDRS}
  ${WIN_RESOURCE_FILE}
  ${CUDA_OBJ_FILES}
  ${cpp_files}
)
QT5_USE_MODULES(${PROJECT_NAME} Widgets OpenGL)
target_link_libraries(${PROJECT_NAME}
#  Qt5::Widgets
#  Qt5::OpenGL
  ${CORE_LIBS}
  ${GLEW_LIBRARIES}
  ${CUDA_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${GLUT_LIBRARIES}
)
