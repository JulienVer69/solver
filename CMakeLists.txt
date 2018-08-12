cmake_minimum_required(VERSION 3.10.2)
#cmakedefine USE_MKL
project(prog)
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect} ")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

#
# Compile.
#
file(GLOB_RECURSE sources  src/*.f90 )
add_executable(linear_solver ${sources} )

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check ")
    add_compile_options(-fdefault-integer-8 -m64 -fopenmp)
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()


find_package(MKL)
if(MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIRS})
    message("-- MKL include ${MKL_INCLUDE_DIRS}") 
    target_link_libraries(linear_solver ${MKL_LIBRARIES} -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -fopenmp)
    message("-- MKL lib ${MKL_LIBRARIES}")

else()
    find_package(BLAS REQUIRED)
    target_link_libraries(linear_solver blas)
    endif()



