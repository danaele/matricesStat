cmake_minimum_required(VERSION 2.8)
 
PROJECT(matricesStat)

# Find SlicerExecutionModel
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

#Find VTK
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

##Find ITK
#find_package(ITK REQUIRED)
#include(${ITK_USE_FILE})

set(matricesStat_SOURCE matricesStat.cxx)

set(matricesStat_HEADER matricesStat.h)

SEMMacroBuildCLI(
  NAME matricesStat
  EXECUTABLE_ONLY
  ADDITIONAL_SRCS ${matricesStat_SOURCE}
  TARGET_LIBRARIES ${VTK_LIBRARIES}
  )

