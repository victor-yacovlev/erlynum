#ifndef NVECTOR_H
#define NVECTOR_H

#include "nblas.h"
#include "narray.h"
#include "ntypes.h"

#include <erl_nif.h>

typedef struct nvector {
    view_params_t   view;
    narray_t        array;
} nvector_t;


ERL_NIF_TERM
erl_nvector_full(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
nvector_to_erl_record(ErlNifEnv * env, nvector_t * vec, const ERL_NIF_TERM bin_ref);

ERL_NIF_TERM
erl_nvector_from_list(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nvector_to_list(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

scalar_value_t
nvector_get(const nvector_t *vec, size_t index);

ERL_NIF_TERM
erl_nvector_get(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nvector_range(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nvector_linspace(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nvector_logspace(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nvector_geomspace(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nvector_copy(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nvector_scale(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nvector_axpy(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]);


#endif // NVECTOR_H
