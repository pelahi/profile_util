cmake_minimum_required(VERSION 3.18)

project(profile_util)

# define cxx 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)
# add cmake macros
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include (PUCMakeUtils)

macro(pu_option optname optdesc status)
   option(PU_${optname} "${optdesc}" "${status}")
endmacro()
pu_option(ENABLE_MPI "Enable mpi" ON)
pu_option(ENABLE_OPENMP "Enable OpenMP" ON)
#pu_option(ENABLE_CUDA "Enable CUDA" OFF)
#pu_option(ENABLE_ROCM "Enable ROCM" OFF)


# check for MPI 
pu_mpi()

# # check for openmp 
pu_openmp()

# #check for cuda 
# pu_cuda()

# # check for rocm 
# pu_rocm()

# set the paths for building
set(PU_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(PU_SOURCE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
include_directories(${PU_INCLUDE_DIRS})
include_directories(${PU_SOURCE_DIRS})
add_subdirectory(${PU_SOURCE_DIRS})
