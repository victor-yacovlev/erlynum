#ifndef ERL_UTIL_H
#define ERL_UTIL_H


#include "ntypes.h"
#include "nvector.h"
#include "nmatrix.h"
#include <erl_nif.h>

ERL_NIF_TERM
make_atom(ErlNifEnv * env, const char *atom);

ERL_NIF_TERM
make_error(ErlNifEnv * env, const char *error);

ERL_NIF_TERM
make_error_with_reason(ErlNifEnv * env, const char *error, const char *reason);

ERL_NIF_TERM
make_dtype(ErlNifEnv * env, const dtype_t dtype);

ERL_NIF_TERM
make_shape(ErlNifEnv * env, size_t shape0, size_t shape1);

ERL_NIF_TERM
make_view_params(ErlNifEnv * env, const view_params_t * params);

ERL_NIF_TERM
make_scalar_value(ErlNifEnv * env, const scalar_element_t scalar, const convert_option_t convert);


_Bool
parse_bool(ErlNifEnv * env, const ERL_NIF_TERM term, _Bool *out, const char* *error);

_Bool
parse_erl_scalar_term(ErlNifEnv * env,
                      const ERL_NIF_TERM term,
                      scalar_element_t *scalar,
                      const create_options_t *options,
                      const char* *error);

_Bool
parse_erl_number(ErlNifEnv * env, const ERL_NIF_TERM term,
                 double * out, _Bool * may_be_single_precision,
                 const char* *error);

_Bool
parse_range_create_options(ErlNifEnv * env, const ERL_NIF_TERM term,
                           range_create_options_t * out,
                           const range_create_options_t * defaults,
                           const char* *error);

_Bool
parse_convert_option(ErlNifEnv * env, const ERL_NIF_TERM term,
                     convert_option_t * out,
                     const char* *error);

_Bool
parse_create_options(ErlNifEnv * env, const ERL_NIF_TERM term,
                     create_options_t * out,
                     const create_options_t * defaults,
                     const char* *error);

_Bool
parse_erl_shape(ErlNifEnv * env, const ERL_NIF_TERM term,
            const size_t arity, size_t shape[arity], const char* *error);

_Bool
parse_view_params(ErlNifEnv * env, const ERL_NIF_TERM term,
                  view_params_t * params, const char* *error);

_Bool
parse_erl_vector(ErlNifEnv * env, const ERL_NIF_TERM term,
                 nvector_t * view, const char* *error,
                 ERL_NIF_TERM *erl_bin_ref);


_Bool
parse_erl_array(ErlNifEnv * env, const ERL_NIF_TERM term,
                 narray_t * arr, const size_t arr_arity,
                 const char* *error,
                 ERL_NIF_TERM *erl_bin_ref);


_Bool
parse_erl_matrix(ErlNifEnv * env, const ERL_NIF_TERM term,
                 nmatrix_t * view, const char* *error,
                 ERL_NIF_TERM *erl_bin_ref);

size_t
min_size(size_t a, size_t b);

size_t
max_size(size_t a, size_t b);


#endif // ERL_UTIL_H
