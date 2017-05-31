#include "erl_util.h"
#include "error_atoms.h"
#include "ntypes.h"
#include "nvector.h"
#include "nmatrix.h"

#include <float.h>
#include <stdbool.h>
#include <string.h>

ERL_NIF_TERM
make_atom(ErlNifEnv * env, const char *atom)
{
    ERL_NIF_TERM ret;
    if (!enif_make_existing_atom(env, atom, &ret, ERL_NIF_LATIN1)) {
        return enif_make_atom(env, atom);
    }
    return ret;
}


ERL_NIF_TERM
make_error(ErlNifEnv * env, const char *error)
{
    ERL_NIF_TERM error_atom = make_atom(env, "error");
    ERL_NIF_TERM message_atom = make_atom(env, error);
    ERL_NIF_TERM tuple = enif_make_tuple2(env, error_atom, message_atom);
    return tuple;
}

ERL_NIF_TERM
make_error_with_reason(ErlNifEnv * env, const char *error, const char *reason)
{
    ERL_NIF_TERM error_atom = make_atom(env, "error");
    ERL_NIF_TERM message_atom = make_atom(env, error);
    ERL_NIF_TERM reason_string = enif_make_string(env, reason, ERL_NIF_LATIN1);
    ERL_NIF_TERM tuple = enif_make_tuple2(env,
                                          error_atom,
                                          enif_make_tuple2(env,
                                                           message_atom,
                                                           reason_string
                                                           )
                                          );
    return tuple;
}


_Bool
parse_erl_scalar_term(ErlNifEnv *env, const ERL_NIF_TERM term,
                      scalar_element_t *scalar,
                      const create_options_t *options,
                      const char **error)
{
    _Bool may_be_single_precision = 0;    
    _Bool force_single_precision = DTSingle==options->dtype || DTComplex==options->dtype || PSingle==options->precision;
    _Bool force_complex = DTComplex==options->dtype || DTDoubleComplex==options->dtype;
    if (enif_is_number(env, term)) {
        double dval = 0;
        if (!parse_erl_number(env, term, &dval, &may_be_single_precision, error)) {
            return false;
        }
        if ( force_single_precision ) {
            if (force_complex) {
                scalar->value.c = CMPLXF((float) dval, 0.0);
                scalar->dtype = DTComplex;
            }
            else {
                scalar->value.s = (float) dval;
                scalar->dtype = DTSingle;
            }

        }
        else {
            if (force_complex) {
                scalar->value.z = CMPLX(dval, 0.0);
                scalar->dtype = DTDoubleComplex;
            }
            else {
                scalar->value.d = dval;
                scalar->dtype = DTDouble;
            }
        }
        return true;
    }
    else if (enif_is_tuple(env, term)) {
        int arity = 0;
        const ERL_NIF_TERM * tuple_items = 0;
        enif_get_tuple(env, term, &arity, &tuple_items);
        if (2!=arity || !enif_is_number(env, tuple_items[0]) || !enif_is_number(env, tuple_items[1])) {
            *error = ERR_ARG_BAD_NUMBER_OR_COMPLEX;
            return false;
        }
        double re = 0.0, im = 0.0;
        _Bool re_single = 0, im_single = 0;
        if (!parse_erl_number(env, tuple_items[0], &re, &re_single, error) ||
                !parse_erl_number(env, tuple_items[1], &im, &im_single, error))
        {
            return false;
        }
        if ( force_single_precision )
        {
            scalar->value.c = CMPLXF((float) re, (float) im);
            scalar->dtype = DTComplex;
        }
        else {
            scalar->value.z = CMPLX(re, im);
            scalar->dtype = DTDoubleComplex;
        }
        return true;
    }
    else {
        *error = ERR_ARG_BAD_NUMBER_OR_COMPLEX;
        return false;
    }
}

_Bool
parse_erl_number(ErlNifEnv *env, const ERL_NIF_TERM term,
                 double *out, _Bool *may_be_single_precision,
                 const char **error)
{
    uint64_t u64 = 0;
    int64_t i64 = 0;
    double dval = 0;
    float fval = 0;
    if (enif_get_uint64(env, term, &u64)) {
        i64 = (int64_t) u64;
        fval = (float) u64;
        dval = (double) u64;
    }
    else if (enif_get_int64(env, term, &i64)) {
        fval = (float) i64;
        dval = (double) i64;
    }
    else if (enif_get_double(env, term, &dval)) {
        fval = (float) fval;
    }
    else {
        *error = ERR_ARG_BAD_NUMBER;
        return false;
    }
    *out = dval;
    *may_be_single_precision = maybe_single_precision(dval);
    return true;
}


_Bool
parse_bool(ErlNifEnv * env, const ERL_NIF_TERM term, _Bool *out, const char* *error)
{
    char value[16]; memset(value, 0, sizeof(value));
    if (!enif_get_atom(env, term, value, sizeof(value), ERL_NIF_LATIN1)) {
        *error = ERR_ARG_BAD_BOOL;
        return false;
    }
    if (0==strncmp("true", value, sizeof(value))) {
        *out = true;
        return true;
    }
    else if (0==strncmp("false", value, sizeof(value))) {
        *out = false;
        return true;
    }
    else {
        *error = ERR_ARG_BAD_BOOL;
        return false;
    }
}

static
_Bool
parse_dtype_atom(const char* atom_value, const size_t max_len, dtype_t *dtype)
{
    if (0==strncmp("auto", atom_value, max_len))
    {
        *dtype = DTAuto;
        return true;
    }
    if (0==strncmp("s", atom_value, max_len) ||
            0==strncmp("single", atom_value, max_len) ||
            0==strncmp("float", atom_value, max_len) ||
            0==strncmp("float32", atom_value, max_len))
    {
        *dtype = DTSingle;
        return true;
    }
    if (0==strncmp("d", atom_value, max_len) ||
            0==strncmp("double", atom_value, max_len) ||
            0==strncmp("float64", atom_value, max_len))
    {
        *dtype = DTDouble;
        return true;
    }
    if (0==strncmp("c", atom_value, max_len) ||
            0==strncmp("complex", atom_value, max_len) ||
            0==strncmp("single_complex", atom_value, max_len) ||
            0==strncmp("float_complex", atom_value, max_len))
    {
        *dtype = DTComplex;
        return true;
    }
    if (0==strncmp("z", atom_value, max_len) ||
            0==strncmp("zcomplex", atom_value, max_len) ||
            0==strncmp("double_complex", atom_value, max_len))
    {
        *dtype = DTDoubleComplex;
        return true;
    }
    return false;
}

static
_Bool
parse_precision_atom(ErlNifEnv * env, const char* atom_value, const size_t max_len, precision_t *prec)
{
    if (0==strncmp("single", atom_value, max_len))
    {
        *prec = PSingle;
        return true;
    }
    if (0==strncmp("double", atom_value, max_len))
    {
        *prec = PDouble;
        return true;
    }
    return false;
}

static _Bool
check_in_coma_separated_string(const char* coma_separated_list, const char* s)
{
    if (!coma_separated_list || !s)
        return false;
    const char * p = strstr(coma_separated_list, s);
    if (!p || '\0'==*p) {
        return false;
    }
    size_t s_len = strlen(s);
    p += s_len;
    return '\0'==*p || ','==*p;
}

static
_Bool
parse_one_create_option(ErlNifEnv * env, const ERL_NIF_TERM term,
                        create_options_t *out,
                        const char* prop_keys_to_ignore,
                        const char* *error)
{
    if (enif_is_tuple(env, term)) {
        int arity = 0;
        const ERL_NIF_TERM * items = 0;
        enif_get_tuple(env, term, &arity, &items);
        if (2!=arity || !enif_is_atom(env, items[0]) || !enif_is_atom(env, items[1])) {
            *error = ERR_ARG_BAD_CREATE_OPTION;
            return false;
        }
        char atom_key[16];      memset(atom_key, 0, sizeof(atom_key));
        char atom_value[16];    memset(atom_value, 0, sizeof(atom_value));
        enif_get_atom(env, items[0], atom_key, sizeof(atom_key), ERL_NIF_LATIN1);
        enif_get_atom(env, items[1], atom_value, sizeof(atom_value), ERL_NIF_LATIN1);
        if (0==strncmp("dtype", atom_key, sizeof(atom_key))) {
            if (!parse_dtype_atom(atom_value, sizeof(atom_value), &out->dtype))
            {
                *error = ERR_ARG_BAD_DTYPE;
                return false;
            }
        }
        else if (0==strncmp("precision", atom_key, sizeof(atom_key))){
            if (!parse_precision_atom(env, atom_value, sizeof(atom_value), &out->precision)) {
                *error = ERR_ARG_BAD_PRECISION;
                return false;
            }
        }
        else if (!check_in_coma_separated_string(prop_keys_to_ignore, atom_key)) {
            *error = ERR_ARG_BAD_CREATE_OPTION;
            return false;
        }
        return true;
    }
    else if (enif_is_atom(env, term)) {
        char atom_value[16];    memset(atom_value, 0, sizeof(atom_value));
        if (!parse_precision_atom(env, atom_value, sizeof(atom_value), &out->precision) &&
                !parse_dtype_atom(atom_value, sizeof(atom_value), &out->dtype)
                )
        {
            char atom_key[16];      memset(atom_key, 0, sizeof(atom_key));
            enif_get_atom(env, term, atom_key, sizeof(atom_key), ERL_NIF_LATIN1);
            if (!check_in_coma_separated_string(prop_keys_to_ignore, atom_key)) {
                *error = ERR_ARG_BAD_CREATE_OPTION;
                return false;
            }
        }
        return true;
    }
    else {
        *error = ERR_ARG_BAD_CREATE_OPTION;
        return false;
    }
}


_Bool
find_boolean_in_prop_list(ErlNifEnv *env,
                          const ERL_NIF_TERM term,
                          const char *key_name,
                          bool *out,
                          const bool defaultValue,
                          const char **error)
{
    *error = 0;
    *out = defaultValue;
    if (!enif_is_list(env, term)) {
        *error = ERR_ARG_BAD_PROPLIST;
        return false;
    }
    ERL_NIF_TERM hd = 0, tl = term;
    char key_atom[64];
    while (enif_get_list_cell(env, tl, &hd, &tl)) {
        memset(key_atom, 0, sizeof(key_atom));
        if (enif_is_atom(env, hd)) {
            enif_get_atom(env, hd, key_atom, sizeof(key_atom), ERL_NIF_LATIN1);
            if (0==strncmp(key_atom, key_name, sizeof(key_atom))) {
                *out = true;
                return true;
            }
        }
        else if (enif_is_tuple(env, hd)) {
            int tuple_arity = 0;
            const ERL_NIF_TERM * tuple_items = 0;
            enif_get_tuple(env, hd, &tuple_arity, &tuple_items);
            if (tuple_arity<2 || !enif_get_atom(env, tuple_items[0], key_atom, sizeof(key_atom), ERL_NIF_LATIN1)) {
                *error = ERR_ARG_BAD_PROPLIST;
                return false;
            }
            if (0==strncmp(key_atom, key_name, sizeof(key_atom))) {
                return parse_bool(env, tuple_items[1], out, error);
            }
        }
    }
    return true;
}

_Bool
parse_create_options(ErlNifEnv *env, const ERL_NIF_TERM term,
                     create_options_t *out,
                     const create_options_t *defaults,
                     const char* prop_keys_to_ignore,
                     const char **error)
{
    if (defaults) {
        memcpy(out, defaults, sizeof(*out));
    }
    else {
        memset(out, 0, sizeof(*out));
    }
    _Bool status = true;
    if (enif_is_list(env, term)) {
        ERL_NIF_TERM hd = 0, tl = term;
        while (enif_get_list_cell(env, tl, &hd, &tl)) {
            status = parse_one_create_option(env, hd, out, prop_keys_to_ignore, error) && status;
        }
    }
    else {
        *error = ERR_ARG_BAD_PROPLIST;
        return false;
    }
    return status;
}

_Bool
parse_erl_shape(ErlNifEnv *env, const ERL_NIF_TERM term,
            const size_t arity, size_t shape[],
            const char **error)
{
    int tuple_arity = 0;
    unsigned list_length = 0;
    size_t provided_arity = 0;
    ERL_NIF_TERM first_term = 0, second_term = 0;
    const ERL_NIF_TERM * tuple_items = 0;
    if (enif_get_tuple(env, term, &tuple_arity, &tuple_items)) {
        if (tuple_arity > 0)
            first_term = tuple_items[0];
        if (tuple_arity > 1)
            second_term = tuple_items[1];
        provided_arity = (size_t) tuple_arity;
    }
    else if (enif_is_list(env, term)) {
        enif_get_list_length(env, term, &list_length);
        provided_arity = (size_t) list_length;
        ERL_NIF_TERM hd = 0, tl = term;
        int i = 0;
        while (enif_get_list_cell(env, tl, &hd, &tl)) {
            if (0==i)
                first_term = hd;
            else if (1==i)
                second_term = hd;
            else
                break;
            ++i;
        }
    }
    else if (enif_is_number(env, term)) {
        provided_arity = 1;
        first_term = term;
        second_term = term;
    }
    else {
        *error = ERR_ARG_BAD_SHAPE;
        return false;
    }
    if (0==provided_arity) {
        *error = ERR_ARG_BAD_SHAPE;
        return false;
    }
    else if (provided_arity < arity) {
        second_term = first_term;
    }

    const ERL_NIF_TERM terms[] = { first_term, second_term };
    for (size_t i=0; i<arity; ++i) {
        unsigned uvalue = 0;
        if (!enif_get_uint(env, terms[i], &uvalue)) {
            *error = ERR_ARG_BAD_SHAPE;
            return false;
        }
        shape[i] = uvalue;
    }
    return true;
}

static _Bool
parse_view_params_one_tuple(ErlNifEnv * env,
                            const ERL_NIF_TERM offset_term,
                            const ERL_NIF_TERM size_term,
                            const ERL_NIF_TERM increment_term,
                            view_params_t * params)
{
    unsigned offset = 0, size = 0;
    int increment = 0;
    if (enif_get_uint(env, offset_term, &offset) &&
            enif_get_uint(env, size_term, &size) &&
            enif_get_int(env, increment_term, &increment))
    {
            params->increment = increment;
            params->offset = offset;
            params->size = size;
            return true;
    }
    return false;
}

_Bool
parse_view_params(ErlNifEnv * env, const ERL_NIF_TERM term,
                  view_params_t *params,
                  const char* *error)
{
    if (enif_is_tuple(env, term)) {
        int tuple_arity = 0;
        const ERL_NIF_TERM * tuple_items = 0;
        enif_get_tuple(env, term, &tuple_arity, &tuple_items);
        char record_name[16];   memset(record_name, 0, sizeof(record_name));
        if (4==tuple_arity && enif_get_atom(env, tuple_items[0], record_name, sizeof(record_name), ERL_NIF_LATIN1))
        {
            if (0==strncmp("view_params", record_name, sizeof(record_name))) {
                return parse_view_params_one_tuple(env, tuple_items[1], tuple_items[2], tuple_items[3], params);
            }

        }
    }
    *error = ERR_REC_VIEW;
    return false;
}


ERL_NIF_TERM
make_dtype(ErlNifEnv *env, const dtype_t dtype)
{
    switch (dtype) {
    case DTSingle:          return make_atom(env, "s");
    case DTDouble:          return make_atom(env, "d");
    case DTComplex:         return make_atom(env, "c");
    case DTDoubleComplex:   return make_atom(env, "z");
    case DTAuto:            return make_atom(env, "auto");
    default:                return make_error(env, "unknown_internal_dtype_value");
    }
}

ERL_NIF_TERM
make_shape(ErlNifEnv *env, size_t shape0, size_t shape1)
{
    return enif_make_tuple2(env,
                            enif_make_uint(env, shape0),
                            enif_make_uint(env, shape1)
                            );
}

ERL_NIF_TERM
make_scalar_value(ErlNifEnv *env, const scalar_element_t scalar, const convert_option_t convert)
{
    ERL_NIF_TERM result = 0;
    double re = 0;
    double im = 0;
    convert_option_t cvt = convert;
    switch (scalar.dtype) {
    case DTSingle:
        re = scalar.value.s;
        cvt = CVTnoconvert==convert ? CVTreal : convert;
        break;
    case DTDouble:
        re = scalar.value.d;
        cvt = CVTnoconvert==convert ? CVTreal : convert;
        break;
    case DTComplex:
        re = crealf(scalar.value.c);
        im = cimagf(scalar.value.c);
        cvt = CVTnoconvert==convert ? CVTcomplex : convert;
        break;
    case DTDoubleComplex:
        re = creal(scalar.value.z);
        im = cimag(scalar.value.z);
        cvt = CVTnoconvert==convert ? CVTcomplex : convert;
        break;
    default:
        break;
    }
    if (CVTinteger==cvt) {
        if (re >= 0) {
            uint64_t uval = (uint64_t) re;
            result = enif_make_uint64(env, uval);
        }
        else {
            int64_t ival = (int64_t) re;
            result = enif_make_int64(env, ival);
        }
    }
    else if (CVTreal==cvt) {
        result = enif_make_double(env, re);
    }
    else if (CVTcomplex==cvt) {
        result = enif_make_tuple2(env,
                                  enif_make_double(env, re),
                                  enif_make_double(env, im)
                                  );
    }
    else {
        result = make_atom(env, "NaN");
    }
    return result;
}

_Bool
parse_erl_vector(ErlNifEnv *env, const ERL_NIF_TERM term,
                 nvector_t *vec, const char **error,
                 ERL_NIF_TERM *erl_bin_ref)
{
    int arity = 0;
    const ERL_NIF_TERM * tuple_items = 0;
    if (!enif_get_tuple(env, term, &arity, &tuple_items) || 5!=arity) {
        *error = ERR_REC_NVECTOR;
        return false;
    }
    char record_name[16];
    if (5!=arity || !enif_get_atom(env, tuple_items[0], record_name, sizeof(record_name), ERL_NIF_LATIN1)) {
        *error = ERR_REC_NVECTOR;
        return false;
    }
    if (0==strncmp("nvector", record_name, sizeof(record_name))) {
        if (!parse_erl_array(env, term, &vec->array, 1, error, erl_bin_ref)) {
            return false;
        }
        if (!parse_view_params(env, tuple_items[1], &vec->view, error)) {
            return false;
        }
        return true;
    }    
    else {
        *error = ERR_REC_NVECTOR;
        return false;
    }
}

_Bool
parse_erl_matrix(ErlNifEnv *env, const ERL_NIF_TERM term,
                 nmatrix_t *mat, const char **error,
                 ERL_NIF_TERM *erl_bin_ref)
{
    int arity = 0;
    const ERL_NIF_TERM * tuple_items = 0;
    if (!enif_get_tuple(env, term, &arity, &tuple_items) || 5!=arity) {
        *error = ERR_REC_NMATRIX;
        return false;
    }
    char record_name[16];
    if (5!=arity || !enif_get_atom(env, tuple_items[0], record_name, sizeof(record_name), ERL_NIF_LATIN1)) {
        *error = ERR_REC_NMATRIX;
        return false;
    }
    if (0==strncmp("nmatrix", record_name, sizeof(record_name))) {
        if (!parse_erl_array(env, term, &mat->array, 2, error, erl_bin_ref)) {
            return false;
        }
        int views_arity = 0;
        const ERL_NIF_TERM * views_tuple = 0;
        if (!enif_get_tuple(env, tuple_items[1], &views_arity, &views_tuple) || 2!=views_arity) {
            *error = ERR_REC_NMATRIX;
            return false;
        }
        if (!parse_view_params(env, views_tuple[0], &mat->view0, error)) {
            return false;
        }
        if (!parse_view_params(env, views_tuple[1], &mat->view1, error)) {
            return false;
        }
        return true;
    }
    else {
        *error = ERR_REC_NMATRIX;
        return false;
    }
}

_Bool
parse_erl_array(ErlNifEnv *env, const ERL_NIF_TERM term,
                narray_t *arr, const size_t arr_arity,
                const char **error,
                ERL_NIF_TERM *erl_bin_ref)
{
    memset(arr, 0, sizeof(*arr));
    int arity = 0;
    const ERL_NIF_TERM * tuple_items = 0;
    enif_get_tuple(env, term, &arity, &tuple_items);

    char dtype_name[16];
    if (!enif_get_atom(env, tuple_items[3], dtype_name, sizeof(dtype_name), ERL_NIF_LATIN1)) {
        *error = ERR_REC_DTYPE;
        return false;
    }
    if (!parse_dtype_atom(dtype_name, sizeof(dtype_name), &arr->dtype) || DTAuto==arr->dtype) {
        *error = ERR_REC_DTYPE;
        return false;
    }
    size_t sizes[2];
    if (!parse_erl_shape(env, tuple_items[2], 2, sizes, error)) {
        return false;
    }
    arr->size0 = sizes[0];
    arr->size1 = sizes[1];

    if (!enif_inspect_binary(env, tuple_items[4], &arr->bin)) {
        *error = ERR_REC_BIN;
        return false;
    }

    if (erl_bin_ref) {
        *erl_bin_ref = tuple_items[4];
    }
    return true;
}



_Bool
parse_convert_option(ErlNifEnv *env, const ERL_NIF_TERM term,
                     convert_option_t *out, const char **error)
{
    char buffer[16];
    memset(buffer, 0, sizeof(buffer));
    if (!enif_get_atom(env, term, buffer, sizeof(buffer), ERL_NIF_LATIN1)) {
        *error = ERR_ARG_BAD_CONVERT;
        return false;
    }
    if (0==strncmp("noconvert", buffer, sizeof(buffer))) {
        *out = CVTnoconvert;
        return true;
    }
    else if (0==strncmp("integer", buffer, sizeof(buffer))) {
        *out = CVTinteger;
        return true;
    }
    else if (0==strncmp("real", buffer, sizeof(buffer))) {
        *out = CVTreal;
        return true;
    }
    else if (0==strncmp("complex", buffer, sizeof(buffer))) {
        *out = CVTcomplex;
        return true;
    }
    else {
        *error = ERR_ARG_BAD_CONVERT;
        return false;
    }
}


ERL_NIF_TERM
make_view_params(ErlNifEnv *env, const view_params_t *params)
{

    return enif_make_tuple4(env,
                            make_atom(env, "view_params"),
                            enif_make_uint(env, params->offset),
                            enif_make_uint(env, params->size),
                            enif_make_int(env, params->increment)
                            );
}

size_t
min_size(size_t a, size_t b)
{
    return a < b ? a : b;
}

size_t
max_size(size_t a, size_t b)
{
    return a > b ? a : b;
}


