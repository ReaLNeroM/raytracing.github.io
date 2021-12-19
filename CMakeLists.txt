#---------------------------------------------------------------------------------------------------
# CMake Build Configuration for the Ray Tracing Weekend Series
#
# See README.md for guidance.
#---------------------------------------------------------------------------------------------------

cmake_minimum_required ( VERSION 3.1.0 )

project ( RTWeekend
  VERSION 3.0.0
  LANGUAGES CXX
)

# Set to c++17
set ( CMAKE_CXX_STANDARD 17 )
set ( CMAKE_CXX_FLAGS "-g" )

# Source
set ( COMMON_ALL
  src/common/rtweekend.h
  src/common/camera.h
  src/common/ray.h
  src/common/vec3.h
)

set ( SOURCE_ORIGINAL
  ${COMMON_ALL}
  src/Original/hittable.h
  src/Original/hittable_list.h
  src/Original/material.h
  src/Original/sphere.h
  src/Original/main.cc
)

set ( SOURCE_BVH
  ${COMMON_ALL}
  src/BVH/hittable.h
  src/BVH/hittable_list.h
  src/BVH/material.h
  src/BVH/sphere.h
  src/BVH/main.cc
)

# Executables
add_executable(render_original      ${SOURCE_ORIGINAL})
add_executable(render_bvh           ${SOURCE_BVH})

include_directories(src/common)
