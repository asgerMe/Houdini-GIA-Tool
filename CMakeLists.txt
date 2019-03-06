# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.6 )
project( My_HDK_Project )
# Locate Houdini's libraries and header files.
# Registers an imported library target named 'Houdini'.
list( APPEND CMAKE_PREFIX_PATH "C:/Program Files/Side Effects Software/Houdini 17.0.459/toolkit/cmake" )

find_package( Houdini REQUIRED )
# Add a library with two source files.
set( library_name FEM_Visco )
add_library( ${library_name} SHARED
    SIM_FEMSolver.cpp
    SIM_FEMSolver.h
    MaterialProperties.h
    GetFEMMatrices.cpp
    GetFEMMatrices.h
    ConstraintMatrix.cpp
    ConstraintMatrix.h
    Placeholder.cpp
    Placeholder.h
    Collision.cpp
    Collision.h
    Solve_System.cpp
    Solve_System.h
    Collision.h
)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.45.0 COMPONENTS shared_ptr.hpp) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    add_executable(progname file1.cxx file2.cxx) 
    target_link_libraries(progname ${Boost_LIBRARIES})
endif()
# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${library_name} Houdini )
# Configure several common target properties, such as its output directory.
houdini_configure_target( ${library_name} )