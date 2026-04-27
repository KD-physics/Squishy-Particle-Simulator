cmake_minimum_required(VERSION 3.18)
project(epd_candidacy CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimisation flags
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# ── pybind11 ──────────────────────────────────────────────────────────────────
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_candidacy_cpp
    src/cpp/bindings.cpp
)

target_include_directories(_candidacy_cpp PRIVATE src/cpp)

install(TARGETS _candidacy_cpp DESTINATION src/simulation)
