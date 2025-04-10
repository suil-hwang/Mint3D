cmake_minimum_required(VERSION 3.12)
project(02_SymFrameView)

# Find Eigen
find_package(Eigen3 3.3 REQUIRED)
# find_package(OpenMP REQUIRED)

file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp ../../src/CubeCover/*.h ../../src/CubeCover/*.cpp ../../src/CubeCover/*.hpp)
# file(GLOB FRONTENDFILES 
#     src/mint3D_hook/*.cpp 
#     src/mint3D_hook/*.h 
#     src/mint3D_hook/*.hpp 
#     src/SymmetricTensor/*.h 
#     src/SymmetricTensor/*.cpp 
#     src/MiNT3D/*.cpp 
#     src/MiNT3D/*.h
#     src/Optimization/*.h
#     src/Optimization/*.cpp)


add_executable(${PROJECT_NAME}  ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR} ../../src ../../lib src/mint3D_hook)
# target_sources(${PROJECT_NAME} PUBLIC ${FRONTENDFILES})

# Conditionally link against Winsock2 library for Windows
if(WIN32)
    target_link_libraries(${PROJECT_NAME} Ws2_32)
endif()


# Use C++17
target_compile_definitions(${PROJECT_NAME} PUBLIC SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}")

# MIP solvers
if(USE_COMISO)
    target_link_libraries(${PROJECT_NAME} igl_copyleft::comiso)
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAS_COMISO)
endif()
if(USE_GUROBI)
    target_link_libraries(${PROJECT_NAME} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY} )
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAS_GUROBI)
endif()


# https://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html 
# here's a good reason to upgrade to C++ 17
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(${PROJECT_NAME} gpgpt polyscope igl::core args nlohmann_json::nlohmann_json) # OpenMP::OpenMP_CXX)

target_link_libraries(${PROJECT_NAME} gpgpt TBB::tbb)

# add_dependencies(${PROJECT_NAME}  ${FRONTENDFILES})
target_compile_definitions(${PROJECT_NAME} PUBLIC -DEIGEN_STACK_ALLOCATION_LIMIT=8388608)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17 )

# Increase stack size, without this integration crashes on larger meshes
# will potentially cause problems if you machine doesn't have sufficient memory.  
# or might need to be increased to process something really big 
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # For Linux, consider using ulimit in the shell instead of setting stack size via linker flags
    message(WARNING "Setting stack size via linker flags is not supported on Linux. Use ulimit in the shell instead.  e.g. ulimit -s 512000")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    # For macOS, set stack size to 512MB (hexadecimal 512MB = 0x20000000)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-stack_size,0x20000000")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    # For Windows (MSVC), set stack size to 512MB
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:536870912")
else()
    message(WARNING "Unsupported platform. Stack size setting not applied.")
endif()