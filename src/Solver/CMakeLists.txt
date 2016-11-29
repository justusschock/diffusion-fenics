cmake_minimum_required(VERSION 3.1)

target_sources(diffusion-fenics
        PUBLIC
        "${CMAKE_CURRENT_LIST_DIR}/convectionDiffusionSolver.h"
        "${CMAKE_CURRENT_LIST_DIR}/currentSolver.h"
        "${CMAKE_CURRENT_LIST_DIR}/poissonSolver.h"
        )
