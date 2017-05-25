#ifndef BACKEND_LOADER_H
#define BACKEND_LOADER_H

#include <stdbool.h>
#include <erl_nif.h>

_Bool
init_backend(const char* priv_dir, const char* backend_name,
             const char* *error, const char* *error_reason);

void*
resolve_blas_function(const char* func_name);



ERL_NIF_TERM
erl_init_backend(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]);

#endif // BACKEND_LOADER_H
