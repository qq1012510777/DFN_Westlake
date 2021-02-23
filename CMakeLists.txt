cmake_minimum_required (VERSION 3.5)

project (DFN_Westlake)

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin) ## bin is the EXECUTABLE_OUTPUT_PATH

set (SRC_LIST ${PROJECT_SOURCE_DIR}/src/DFN_main.cpp)


#########below use a LIB: OPENMP
set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

FIND_PACKAGE( OpenMP REQUIRED) #find omp
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS"${CMAKE_EXE_LINKER_FLAGS}${OpenMP_EXE_LINKER_FLAGS}")
endif()
#########OPENMP

################## head fild
# find head files of local dfn lib
include_directories (${PROJECT_SOURCE_DIR}/include/DFN_H)

# find head files of Fade_2D
include_directories ($ENV{HOME}/external_lib/fadeRelease_v1.84/include_fade2d)

# find head files of Eigen/Dense
include_directories ($ENV{HOME}/external_lib)

# find head files of gsl lib
include_directories ($ENV{HOME}/external_lib/gsl-2.1)
################## head fild

##################pkg

# find pkg of FADE_2D
find_library(FADE2D_LIB fade2d HINTS $ENV{HOME}/external_lib/fadeRelease_v1.84/lib_ubuntu18.04_x86_64)

#find pkg of GSL_CBLAS
find_library(GSL_CBLAS_LIB gslcblas HINTS $ENV{HOME}/external_lib/cblas/.libs)

#find pkg of GSL
find_library(GSL_LIB gsl HINTS $ENV{HOME}/external_lib/gsl-2.1)

#find pkg of LAPACK
find_library(LAPACK_LIB lapack HINTS $ENV{HOME}/external_lib/lapack-3.5.0/lib)

#find pkg of BLAS
find_library(BLAS_LIB blas HINTS $ENV{HOME}/external_lib/lapack-3.5.0/lib )


######################
add_executable (DFN_main ${SRC_LIST})

target_link_libraries (DFN_main ${FADE2D_LIB})

target_link_libraries (DFN_main ${GSL_LIB})

target_link_libraries (DFN_main ${GSL_CBLAS_LIB})

target_link_libraries (DFN_main ${LAPACK_LIB} "gfortran")

target_link_libraries (DFN_main ${BLAS_LIB} "gfortran")
