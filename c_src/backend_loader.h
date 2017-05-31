#ifndef BACKEND_LOADER_H
#define BACKEND_LOADER_H

#include <stdbool.h>
#include <erl_nif.h>

typedef enum {
    DrvNone         = 0x00,
    DrvAuto         = 0xFF,
    DrvIntelMKL     = 0x01,
    DrvATLAS        = 0x02,
    DrvNetlibBLAS   = 0x03,
} blas_backend_driver_t;

_Bool
init_backend(const char* priv_dir, const char* backend_name,
             const char* *error, const char* *error_reason);

void*
resolve_blas_function(const char* func_name);

blas_backend_driver_t
backend_in_use();


ERL_NIF_TERM
erl_init_backend(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_backend_in_use(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]);

#endif // BACKEND_LOADER_H
