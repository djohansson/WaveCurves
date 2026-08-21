cmake_minimum_required(VERSION 3.16)
project(wavecurves C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# --- SDL3 -------------------------------------------------------------------
# Prefer an installed SDL3; otherwise fetch and build the one from
# https://github.com/libsdl-org/SDL. Set WAVECURVES_FETCH_SDL=OFF to require a
# system SDL3 instead.
option(WAVECURVES_FETCH_SDL "Download SDL3 if it is not already installed" ON)

find_package(SDL3 QUIET CONFIG)

if(NOT SDL3_FOUND AND WAVECURVES_FETCH_SDL)
    include(FetchContent)
    set(SDL_SHARED   ON  CACHE BOOL "" FORCE)
    set(SDL_STATIC   OFF CACHE BOOL "" FORCE)
    set(SDL_TEST_LIBRARY OFF CACHE BOOL "" FORCE)
    set(SDL_EXAMPLES OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(SDL3
        GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
        GIT_TAG        release-3.2.24
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(SDL3)
elseif(NOT SDL3_FOUND)
    message(FATAL_ERROR "SDL3 not found and WAVECURVES_FETCH_SDL is OFF")
endif()

# --- target -----------------------------------------------------------------
add_executable(wavecurves
    src/main.c
    src/wc_base.c
    src/wc_curve.c
    src/wc_render.c
    src/wc_png.c
)

target_include_directories(wavecurves PRIVATE src)
target_link_libraries(wavecurves PRIVATE SDL3::SDL3)

if(MSVC)
    target_compile_options(wavecurves PRIVATE /W4 /fp:fast)
    target_compile_definitions(wavecurves PRIVATE _CRT_SECURE_NO_WARNINGS)
else()
    target_compile_options(wavecurves PRIVATE -Wall -Wextra -ffast-math)
    target_link_libraries(wavecurves PRIVATE m)
endif()

# Keep the SDL3 shared library next to the executable on Windows.
if(WIN32 AND TARGET SDL3-shared)
    add_custom_command(TARGET wavecurves POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:SDL3-shared> $<TARGET_FILE_DIR:wavecurves>)
endif()
