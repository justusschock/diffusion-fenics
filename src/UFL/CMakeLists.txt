cmake_minimum_required(VERSION 3.1)

target_sources(diffusion-fenics
        PUBLIC
        "${CMAKE_CURRENT_LIST_DIR}/convectionDiffusion3D.h"
        "${CMAKE_CURRENT_LIST_DIR}/current3D.h"
        "${CMAKE_CURRENT_LIST_DIR}/poisson3D.h"
        )

add_custom_command(
        OUTPUT
            ${CMAKE_CURRENT_LIST_DIR}/convectionDiffusion3D.h
        DEPENDS
            ${CMAKE_CURRENT_LIST_DIR}/convectionDiffusion3D.ufl
        COMMAND
            ffc -l dolfin convectionDiffusion3D.ufl
        VERBATIM
)

add_custom_command(
        OUTPUT
            ${CMAKE_CURRENT_LIST_DIR}/current3D.h
        DEPENDS
            ${CMAKE_CURRENT_LIST_DIR}/current3D.ufl
        COMMAND
            ffc -l dolfin current3D.ufl
        VERBATIM
)

add_custom_command(
        OUTPUT
            ${CMAKE_CURRENT_LIST_DIR}/poisson3D.h
        DEPENDS
            ${CMAKE_CURRENT_LIST_DIR}/poisson3D.ufl
        COMMAND
            ffc -l dolfin poisson3D.ufl
        VERBATIM
)