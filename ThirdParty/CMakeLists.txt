add_subdirectory(polyscope)
add_subdirectory(mma)

file(GLOB IMPLOT_HEADERS "implot/*.h")
file(GLOB IMPLOT_SOURCES "implot/*.cpp")
add_library(implot ${IMPLOT_HEADERS} ${IMPLOT_SOURCES})
target_link_libraries(implot imgui)
