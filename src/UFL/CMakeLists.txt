cmake_minimum_required(VERSION 3.1)

file(GLOB UFL RELATIVE ${CMAKE_CURRENT_LIST_DIR} "*.ufl")
#target_sources(diffusion-fenics
#        PUBLIC
#        "${CMAKE_CURRENT_LIST_DIR}/convectionDiffusion3D.h"
#        "${CMAKE_CURRENT_LIST_DIR}/current3D.h"
#        "${CMAKE_CURRENT_LIST_DIR}/poisson3D.h"
#        )
foreach(ufl_loop ${UFL})
    get_filename_component(ufl_to_h ${ufl_loop} NAME_WE)    
    list(APPEND h_list ${ufl_to_h}.h)
add_custom_command(
        OUTPUT
            ${CMAKE_CURRENT_LIST_DIR}/${ufl_to_h}.h
        DEPENDS
            ${CMAKE_CURRENT_LIST_DIR}/${ufl_loop}
        COMMAND
            ffc -O -fquadrature_degree=4 -l dolfin ${ufl_loop}
        WORKING_DIRECTORY
             ${CMAKE_CURRENT_LIST_DIR}
        VERBATIM
)
endforeach(ufl_loop)
#MESSAGE(${UFL})
#MESSAGE(${h_list})

add_custom_target(generate_ufl DEPENDS ${h_list})
