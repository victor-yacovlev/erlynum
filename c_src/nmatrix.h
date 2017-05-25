#ifndef NMATRIX_H
#define NMATRIX_H
#include "nblas.h"
#include "ntypes.h"
#include "nvector.h"
#include <erl_nif.h>

typedef struct nmatrix {
    view_params_t   view0;
    view_params_t   view1;
    narray_t        array;
} nmatrix_t;

ERL_NIF_TERM
erl_nmatrix_full(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nmatrix_to_list(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nmatrix_from_list(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nmatrix_eye(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nmatrix_row(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nmatrix_col(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
erl_nmatrix_diag(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

scalar_value_t
nmatrix_get(const nmatrix_t *mat, size_t row_index, size_t col_index);

ERL_NIF_TERM
erl_nmatrix_get(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]);

ERL_NIF_TERM
nmatrix_to_erl_record(ErlNifEnv * env, nmatrix_t * mat, const ERL_NIF_TERM bin_ref);

#endif // NMATRIX_H
