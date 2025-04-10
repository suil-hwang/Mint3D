cmake_minimum_required(VERSION 3.12)
project(00_basic_test)

# Find Eigen
find_package(Eigen3 3.3 REQUIRED)

file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp)	  

add_executable(${PROJECT_NAME} ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR} ../../src)
# add_subdirectory(libs/TinyAD)\

# find_package(OpenMP)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
# endif()

# if(OpenMP_CXX_FOUND OR OPENMP_FOUND)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

# # Link against the OpenMP libraries
# if(APPLE AND EXISTS /usr/local/opt/libomp/lib/libomp.dylib)
#     set(OpenMP_LIBRARIES "/usr/local/opt/libomp/lib/libomp.dylib")
#     set(OpenMP_INCLUDES "/usr/local/opt/libomp/include")

#     include_directories("${OpenMP_INCLUDES}")
#     link_directories("/usr/local/opt/libomp/lib")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")
# endif()

# # Then link these to your target (e.g., replace YourTargetName with your actual target)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(${PROJECT_NAME} gpgpt polyscope OpenMP::OpenMP_CXX)

#     # target_link_libraries(${PROJECT_NAME}  PUBLIC OpenMP::OpenMP_CXX)
# elseif(APPLE)
#     target_link_libraries(${PROJECT_NAME} gpgpt polyscope "${OpenMP_LIBRARIES}")
# endif()

# Use C++17
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(${PROJECT_NAME} gpgpt polyscope) # OpenMP::OpenMP_CXX)
