add_library(ApplicationUtilLib INTERFACE)
target_include_directories(ApplicationUtilLib INTERFACE .)
target_link_libraries(ApplicationUtilLib INTERFACE
  CLI11::CLI11
)

add_executable(parameterize_aligned
parameterize_aligned.cpp
)
target_link_libraries(parameterize_aligned PRIVATE
  PennerLib
  ApplicationUtilLib
)

add_executable(add_trivial_quad_parameterization
add_trivial_quad_parameterization.cpp
)
target_link_libraries(add_trivial_quad_parameterization PRIVATE
  PennerLib
  ApplicationUtilLib
)

add_executable(generate_components
generate_components.cpp
)
target_link_libraries(generate_components PRIVATE
  PennerLib
  ApplicationUtilLib
)

# only build visualization executables if visualization enabled
if (NOT USE_MULTIPRECISION)
if (ENABLE_VISUALIZATION)
    add_executable(view_quad_mesh
    view_quad_mesh.cpp
    )
    target_link_libraries(view_quad_mesh PRIVATE
      PennerLib
      ApplicationUtilLib
    )

    add_executable(view_seamless_uv
    view_seamless_uv.cpp
    )
    target_link_libraries(view_seamless_uv PRIVATE
      PennerLib
      ApplicationUtilLib
    )
endif()
endif()
