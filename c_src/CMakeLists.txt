# This project file is not required to build the module.
# Is is used only by IDEs while in development.
# Use `rebar compile` or `rebar3 compile` from parent directory.

cmake_minimum_required(VERSION 3.0)
project(erlynum)

execute_process(
    COMMAND         erl -noshell -s init stop -eval "io:format(\"~s/erts-~s/include/\", [code:root_dir(), erlang:system_info(version)])."
    OUTPUT_VARIABLE ERTS_INCLUDE_DIR
)

execute_process(
    COMMAND         erl -noshell -s init stop -eval "io:format(\"~s\", [code:lib_dir(erl_interface, include)])."
    OUTPUT_VARIABLE ERL_INTERFACE_INCLUDE_DIR
)

execute_process(
    COMMAND         erl -noshell -s init stop -eval "io:format(\"~s\", [code:lib_dir(erl_interface, lib)])."
    OUTPUT_VARIABLE ERL_INTERFACE_LIB_DIR
)

message(STATUS "ERTS_INCLUDE_DIR: ${ERTS_INCLUDE_DIR}")
message(STATUS "ERL_INTERFACE_INCLUDE_DIR: ${ERL_INTERFACE_INCLUDE_DIR}")
message(STATUS "ERL_INTERFACE_LIB_DIR: ${ERL_INTERFACE_LIB_DIR}")

include_directories(${ERTS_INCLUDE_DIR})
include_directories(${ERL_INTERFACE_INCLUDE_DIR})


set(COMMON_FLAGS "-finline-functions -Werror -Wall -Wmaybe-uninitialized -Wreturn-type -Wnonnull -Wstrict-aliasing")
set(CMAKE_C_FLAGS "-std=c99 -O3 -D_GNU_SOURCE -Wmissing-prototypes ${COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS "-std=c++11 -D_GNU_SOURCE ${COMMON_FLAGS}")

set(SOURCES
    erlynum_drv.c
    backend_loader.c
    backend_loader_posix.c
    backend_loader_win32.c
    erl_util.c
    ntypes.c
    narray.c
    nvector.c
    nmatrix.c
    nblas.c
    blas_fallback.c
)

set(DRIVER_OUTPUT "${PROJECT_NAME}_drv")

add_library(${DRIVER_OUTPUT} SHARED ${SOURCES})
target_link_libraries(${DRIVER_OUTPUT} "-L${ERL_INTERFACE_LIB_DIR} -lerl_interface -lei -lm")
target_compile_options(${DRIVER_OUTPUT} PRIVATE ${GEN_FLAGS})
set_property(TARGET ${DRIVER_OUTPUT} PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET ${DRIVER_OUTPUT} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../priv/")
set_property(TARGET ${DRIVER_OUTPUT} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../priv/")
set_property(TARGET ${DRIVER_OUTPUT} PROPERTY PREFIX "")

