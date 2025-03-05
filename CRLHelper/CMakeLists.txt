set(DEPS tbb tbbmalloc tbbmalloc_proxy mkl_intel_lp64 mkl_sequential mkl_core igl::core igl_restricted::triangle imgui polyscope glfw glad CGAL cholmod umfpack klu nlopt implot mma::mma mma::gcmma)

file(GLOB_RECURSE HEADERS "*.h" "*.hpp")
file(GLOB_RECURSE SOURCES "*.cpp")

add_library(CRLHelper ${HEADERS} ${SOURCES})
target_link_libraries(CRLHelper ${DEPS})
