cmake_minimum_required(VERSION 3.1)

#Pre-build regeneration of Header-Files out of UFL-Files to ensure correctness and framwork compatibility
#add_custom_target(PRE_BUILD_FFC_POISSON ALL ffc -l dolfin ${CMAKE_CURRENT_LIST_DIR}/*.ufl)

target_sources(diffusion-fenics
        PRIVATE #${PRE_BUILD_FFC_POISSON}
        "${CMAKE_CURRENT_LIST_DIR}/poisson1D.h"
        "${CMAKE_CURRENT_LIST_DIR}/poisson2D.h"
        "${CMAKE_CURRENT_LIST_DIR}/poisson3D.h"
        "${CMAKE_CURRENT_LIST_DIR}/poissonSolver.h")
