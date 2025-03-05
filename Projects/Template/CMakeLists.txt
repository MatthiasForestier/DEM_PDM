add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
set(DEPS Half tbb tbbmalloc tbbmalloc_proxy mkl_intel_lp64 mkl_sequential mkl_core igl::core igl_restricted::triangle imgui CGAL suitesparseconfig cholmod metis CRLHelper OpenGL)

file(GLOB_RECURSE HEADERS "include/*.h" "include/*.hpp" "../../CRLHelper/*.h")
file(GLOB_RECURSE SOURCES "src/*.cpp" "../../CRLHelper/*.cpp")

add_executable(Template ${HEADERS} ${SOURCES})
target_link_libraries(Template ${DEPS})
