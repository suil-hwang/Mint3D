cmake_minimum_required(VERSION 3.12)
project(01_Mint3D)

# Find Eigen
find_package(Eigen3 3.3 REQUIRED)
# find_package(OpenMP REQUIRED)

file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp ../../src/CubeCover/*.h ../../src/CubeCover/*.cpp ../../src/CubeCover/*.hpp)
file(GLOB FRONTENDFILES 
    src/mint3D_hook/*.cpp 
    src/mint3D_hook/*.h 
    src/mint3D_hook/*.hpp 
    src/SymmetricTensor/*.h 
    src/SymmetricTensor/*.cpp 
    src/Optimization/*.h
    src/Optimization/*.cpp
    src/MiNT3D/*.cpp 
    src/MiNT3D/*.h)


add_executable(${PROJECT_NAME}  ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR} ../../src ../../lib src/mint3D_hook src/MiNT3D)
target_sources(${PROJECT_NAME} PUBLIC ${FRONTENDFILES})

# Conditionally link against Winsock2 library for Windows
if(WIN32)
    target_link_libraries(${PROJECT_NAME} Ws2_32)
endif()


# Use C++17
target_compile_definitions(${PROJECT_NAME} PUBLIC SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}")


# https://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html 
# here's a good reason to upgrade to C++ 17
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(${PROJECT_NAME} gpgpt polyscope igl::core args) #  nasoq) # OpenMP::OpenMP_CXX)

target_link_libraries(${PROJECT_NAME} gpgpt TBB::tbb)

# add_dependencies(${PROJECT_NAME}  ${FRONTENDFILES})
target_compile_definitions(${PROJECT_NAME} PUBLIC -DEIGEN_STACK_ALLOCATION_LIMIT=8388608)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17 )

# Add AddressSanitizer flags only in debug mode
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${PROJECT_NAME} PRIVATE -fsanitize=address -g)
        target_link_options(${PROJECT_NAME} PRIVATE -fsanitize=address)
    endif()
endif()
