
cmake_minimum_required(VERSION 3.1)

set(PROJECT_NAME diffusion-fenics)
project(${PROJECT_NAME})

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set CMake behavior
cmake_policy(SET CMP0004 OLD)

# Get DOLFIN configuration data (DOLFINConfig.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

#Get VTK configuration data
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
            "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Compiler flags
set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")


# Include directories
include_directories(${DOLFIN_INCLUDE_DIRS})
include_directories(SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
include_directories(SYSTEM ${CMAKE_CURRENT_LIST_DIR}/UFL)
include_directories(${CMAKE_CURRENT_LIST_DIR})

set(SOURCE_FILES main.cpp Solver/rungeKuttaDiffusionEquation.h)# pdeSetupClasses.h)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${build_dir}/bin)

# Executable
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
add_dependencies(${PROJECT_NAME} generate_ufl)
#add_subdirectory(Solver)
add_subdirectory(UFL)


# Target libraries
target_compile_options(${PROJECT_NAME} PUBLIC ${_flags_all})
target_link_libraries(${PROJECT_NAME} ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES})
