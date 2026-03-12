cmake_minimum_required(VERSION 3.16)
project(fin_flutter_core VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Dependencies ──────────────────────────────────────────────────────────────
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# ── Header-only library target ────────────────────────────────────────────────
add_library(fin_flutter_core INTERFACE)
target_include_directories(fin_flutter_core INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(fin_flutter_core INTERFACE Eigen3::Eigen)
target_compile_options(fin_flutter_core INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -O2>
)

# ── Tests ─────────────────────────────────────────────────────────────────────
enable_testing()
add_subdirectory(tests)
