include(FetchContent)
FetchContent_Declare(
  polyscope
  GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
  GIT_TAG 2cd38b7496a1e8ff79f615920de2980295bedef7
)

FetchContent_GetProperties(polyscope)
if(NOT polyscope_POPULATED)
  # Fetch the content using previously declared details
  FetchContent_Populate(polyscope)
  message(STATUS "polyscope_SOURCE_DIR: ${polyscope_SOURCE_DIR}")
  message(STATUS "polyscope_BINARY_DIR: ${polyscope_BINARY_DIR}")
  add_subdirectory(${polyscope_SOURCE_DIR} ${polyscope_BINARY_DIR})
endif()
FetchContent_MakeAvailable(polyscope)


# for now just manually include implot in this gpgpt

# # Set the IMGUI_INCLUDE_DIR
# set(IMGUI_INCLUDE_DIR ${polyscope_SOURCE_DIR}/deps/imgui)

# # Fetch implot
# FetchContent_Declare(
#     implot
#     GIT_REPOSITORY https://github.com/epezent/implot.git
#     GIT_TAG v0.16
# )
# FetchContent_GetProperties(implot)
# if(NOT implot_POPULATED)
#   FetchContent_Populate(implot)
#   add_subdirectory(${implot_SOURCE_DIR} ${implot_BINARY_DIR})
# endif()

# # Here, adjust the include directories for implot if needed
# target_include_directories(implot PUBLIC ${IMGUI_INCLUDE_DIR})


