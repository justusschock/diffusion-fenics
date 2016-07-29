cmake_minimum_required(VERSION 3.1)

#Pre-build regeneration of Header-Files out of UFL-Files to ensure correctness and framwork compatibility
<<<<<<< HEAD
add_custom_target(PRE_BUILD_FFC_DIFFUSION ALL ffc -l dolfin ${CMAKE_CURRENT_LIST_DIR}/*.ufl)

target_sources(diffusion-fenics
        PRIVATE ${PRE_BUILD_FFC_DIFFUSION}
        "${CMAKE_CURRENT_LIST_DIR}/convectionDiffusion1D.h"
        "${CMAKE_CURRENT_LIST_DIR}/convectionDiffusion2D.h"
        "${CMAKE_CURRENT_LIST_DIR}/convectionDiffusion3D.h"
        "${CMAKE_CURRENT_LIST_DIR}/convectionDiffusionSolver.h")
=======
#add_custom_target(PRE_BUILD_FFC_DIFFUSION ALL ffc -l dolfin ${CMAKE_CURRENT_LIST_DIR}/*.ufl)

target_sources(diffusion-fenics
        PRIVATE #${PRE_BUILD_FFC_DIFFUSION}
        "${CMAKE_CURRENT_LIST_DIR}/convectionDiffusion1D.h"
        "${CMAKE_CURRENT_LIST_DIR}/convectionDiffusion2D.h"
        "${CMAKE_CURRENT_LIST_DIR}/convectionDiffusion3D.h"
        "${CMAKE_CURRENT_LIST_DIR}/convectionDiffusionSolver.h")
>>>>>>> remotes/origin/fenics_1.6.x
