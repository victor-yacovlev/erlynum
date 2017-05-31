#include "erl_util.h"
#include "error_atoms.h"
#include "ntypes.h"
#include "nvector.h"

#include <erl_nif.h>
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

ERL_NIF_TERM
erl_nvector_full(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (3!=argc) {
        return enif_make_badarg(env);
    }    

    const char * error = 0;
    create_options_t options;
    if (!parse_create_options(env, argv[2], &options, 0, 0, &error)) {
        return make_error(env, error);
    }

    nvector_t result;   memset(&result, 0, sizeof(result));
    scalar_element_t fill_value;

    if (!parse_erl_scalar_term(env, argv[1], &fill_value, &options, &error)) {
        return make_error(env, error);
    }
    result.array.dtype = fill_value.dtype;
    result.array.size0 = 1;

    if (!parse_erl_shape(env, argv[0], 1, &result.array.size1, &error)) {
        return make_error(env, error);
    }
    result.view.offset = 0;
    result.view.size = result.array.size1;
    result.view.increment = 1;

    const size_t item_size = scalar_size(result.array.dtype);
    const size_t bin_size = item_size * result.array.size1;
    if (!enif_alloc_binary(bin_size, &result.array.bin)) {
        return make_error(env, ERR_GEN_MEMORY_ALLOCATION);
    }
    void * p = result.array.bin.data;
    for (size_t i=0; i<result.array.size1; ++i) {
        memcpy(p, &fill_value.value, item_size);
        p += item_size;
    }

    return nvector_to_erl_record(env, &result, 0);
}


ERL_NIF_TERM
nvector_to_erl_record(ErlNifEnv *env, nvector_t *vec, const ERL_NIF_TERM bin_ref)
{
    const ERL_NIF_TERM header   = make_atom(env, "nvector");
    const ERL_NIF_TERM view     = make_view_params(env, &vec->view);
    const ERL_NIF_TERM shape    = make_shape(env, vec->array.size0, vec->array.size1);
    const ERL_NIF_TERM dtype    = make_dtype(env, vec->array.dtype);
    const ERL_NIF_TERM data     = bin_ref ? bin_ref : enif_make_binary(env, &vec->array.bin);
    return enif_make_tuple5(env, header, view, shape, dtype, data);
}

ERL_NIF_TERM
erl_nvector_from_list(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (2!=argc) {
        return enif_make_badarg(env);
    }

    const char * error = 0;
    create_options_t options;
    if (!parse_create_options(env, argv[1], &options, 0, 0, &error)) {
        return make_error(env, error);
    }

    if (!enif_is_list(env, argv[0])) {
        return make_error(env, ERR_ARG_BAD_LIST);
    }

    ERL_NIF_TERM hd = 0, tl = argv[0];
    unsigned elements_count = 0;
    enif_get_list_length(env, argv[0], &elements_count);
    scalar_element_t * elements = enif_alloc(elements_count*sizeof(*elements));
    scalar_element_t * current = 0;
    unsigned element_index = 0;
    _Bool error_occured = 0;
    while (enif_get_list_cell(env, tl, &hd, &tl)) {
        current = elements + element_index*sizeof(*current);
        if (!parse_erl_scalar_term(env, hd, current, &options, &error)) {
            error_occured = 1;
            break;
        }
        element_index++;
    }
    if (error_occured) {
        enif_free(elements);
        return make_error(env, error);
    }

    nvector_t result;
    memset(&result, 0, sizeof(result));
    result.array.size0 = 1;
    result.array.size1 = elements_count;
    result.array.dtype = common_dtype_for_elements(elements_count, elements);
    result.view.offset = 0;
    result.view.size = elements_count;
    result.view.increment = 1;

    const size_t item_size = scalar_size(result.array.dtype);
    const size_t bin_size = item_size * result.array.size1;
    if (!enif_alloc_binary(bin_size, &result.array.bin)) {
        return make_error(env, ERR_GEN_MEMORY_ALLOCATION);
    }
    void * p = result.array.bin.data;
    for (size_t i=0; i<result.array.size1; ++i) {
        current = elements + i*sizeof(*current);
        scalar_value_t fill_value = convert_type(current->value, current->dtype, result.array.dtype);
        memcpy(p, &fill_value, item_size);
        p += item_size;
    }

    return nvector_to_erl_record(env, &result, 0);
}

ERL_NIF_TERM
erl_nvector_to_list(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (2!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nvector_t vec;
    if (!parse_erl_vector(env, argv[0], &vec, &error, NULL)) {
        return make_error(env, error);
    }
    convert_option_t convert;
    if (!parse_convert_option(env, argv[1], &convert, &error)) {
        return make_error(env, error);
    }
    const size_t items_count = vec.view.size;
    const size_t items_start = vec.view.offset;
    const int items_increment = vec.view.increment;
    const size_t item_size = scalar_size(vec.array.dtype);
    ERL_NIF_TERM * list_items = enif_alloc(items_count * sizeof(ERL_NIF_TERM));
    for (size_t i=0; i<items_count; ++i) {
        const size_t real_index = items_start + i * items_increment;
        void * src = vec.array.bin.data + item_size * real_index;
        scalar_element_t element;
        element.dtype = vec.array.dtype;
        memcpy(&element.value, src, item_size);
        list_items[i] = make_scalar_value(env, element, convert);
    }
    ERL_NIF_TERM result = enif_make_list_from_array(env, list_items, items_count);
    enif_free(list_items);
    return result;
}

scalar_value_t
nvector_get(const nvector_t *vec, size_t index)
{
    scalar_value_t element;
    const size_t item_size = scalar_size(vec->array.dtype);
    const size_t real_index = vec->view.offset + index * vec->view.increment;
    void * src = vec->array.bin.data + item_size * real_index;
    memcpy(&element, src, item_size);
    return element;
}

ERL_NIF_TERM
erl_nvector_get(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (3!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nvector_t vec;
    if (!parse_erl_vector(env, argv[0], &vec, &error, NULL)) {
        return make_error(env, error);
    }
    unsigned index = 0;
    if (!enif_get_uint(env, argv[1], &index)) {
        return make_error(env, ERR_ARG_BAD_INDEX);
    }
    convert_option_t convert;
    if (!parse_convert_option(env, argv[2], &convert, &error)) {
        return make_error(env, error);
    }
    const size_t items_count = vec.view.size;
    if (index >= items_count) {
        return make_error(env, ERR_ARG_INDEX_OUT_OF_RANGE);
    }
    scalar_element_t element;
    element.dtype = vec.array.dtype;
    element.value = nvector_get(&vec, index);
    ERL_NIF_TERM result = make_scalar_value(env, element, convert);

    return result;
}





ERL_NIF_TERM
erl_nvector_range(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (4!=argc) {
        return enif_make_badarg(env);
    }
    range_create_options_t create_options;
    const char* error = 0;
    if (!parse_create_options(env, argv[3], &create_options.create_options, 0, "endpoint", &error)) {
        return make_error(env, error);
    }
    if (!find_boolean_in_prop_list(env, argv[3], "endpoint", &create_options.endpoint, false, &error)) {
        return make_error(env, error);
    }
    scalar_element_t start;
    scalar_element_t stop;
    scalar_element_t step;
    if (!parse_erl_scalar_term(env, argv[0], &start, &create_options.create_options, &error) ||
            !parse_erl_scalar_term(env, argv[1], &stop, &create_options.create_options, &error) ||
            !parse_erl_scalar_term(env, argv[2], &step, &create_options.create_options, &error)
            )
    {
        return make_error(env, error);
    }
    nvector_t result;
    memset(&result, 0, sizeof(result));
    result.array.dtype = DTAuto==create_options.create_options.dtype
            ? common_dtype(common_dtype(start.dtype, stop.dtype), step.dtype)
            : create_options.create_options.dtype;
    if (!narray_allocate_range(&result.array, &start, &stop, &step, create_options.endpoint, &error)) {
        return make_error(env, error);
    }
    result.view.offset = 0;
    result.view.increment = 1;
    result.view.size = result.array.size0 * result.array.size1;
    return nvector_to_erl_record(env, &result, 0);
}

ERL_NIF_TERM
erl_nvector_linspace(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (4!=argc) {
        return enif_make_badarg(env);
    }
    range_create_options_t create_options;

    const char* error = 0;
    if (!parse_create_options(env, argv[3], &create_options.create_options, 0, "endpoint", &error)) {
        return make_error(env, error);
    }
    if (!find_boolean_in_prop_list(env, argv[3], "endpoint", &create_options.endpoint, true, &error)) {
        return make_error(env, error);
    }
    scalar_element_t start;
    scalar_element_t stop;
    if (!parse_erl_scalar_term(env, argv[0], &start, &create_options.create_options, &error) ||
            !parse_erl_scalar_term(env, argv[1], &stop, &create_options.create_options, &error)
            )
    {
        return make_error(env, error);
    }
    unsigned count = 0;
    if (!enif_get_uint(env, argv[2], &count)) {
        return make_error(env, ERR_ARG_BAD_COUNT);
    }
    nvector_t result;
    memset(&result, 0, sizeof(result));

    if (0==count) {
        result.array.dtype = DTAuto==create_options.create_options.dtype
                ? common_dtype(start.dtype, stop.dtype)
                : create_options.create_options.dtype;
        result.array.size0 = 1;
        result.array.size1 = 0;
        result.view.offset = 0;
        result.view.size = 0;
        result.view.increment = 1;
        if (!enif_alloc_binary(0, &result.array.bin)) {
            return make_error(env, ERR_GEN_MEMORY_ALLOCATION);
        }
        return nvector_to_erl_record(env, &result, 0);
    }

    if (DTAuto==create_options.create_options.dtype && count > 1) {
        double step_re = double_value(&stop, false) - double_value(&start, false);
        double step_im = double_value(&stop, true) - double_value(&start, true);
        step_re /= (count-1);
        step_im /= (count-1);
        const _Bool single_re = maybe_single_precision(step_re);
        const _Bool single_im = maybe_single_precision(step_im);

        scalar_element_t step;
        if (0==step_im) {
            step.dtype = single_re ? DTSingle : DTDouble;
        }
        else {
            step.dtype = single_re && single_im ? DTComplex : DTDoubleComplex;
        }
        result.array.dtype = common_dtype(common_dtype(start.dtype, stop.dtype), step.dtype);
    }
    else if (DTAuto==create_options.create_options.dtype && 1==count) {
        result.array.dtype = stop.dtype;
    }
    else {
        result.array.dtype = create_options.create_options.dtype;
    }

    result.array.size0 = 1;
    result.array.size1 = count;
    result.view.offset = 0;
    result.view.size = count;
    result.view.increment = 1;

    if (!narray_allocate_linspace(&result.array, &start, &stop, create_options.endpoint, &error)) {
        return make_error(env, error);
    }
    return nvector_to_erl_record(env, &result, 0);
}

ERL_NIF_TERM
erl_nvector_logspace(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (5!=argc) {
        return enif_make_badarg(env);
    }
    range_create_options_t create_options;

    const char* error = 0;
    if (!parse_create_options(env, argv[4], &create_options.create_options, 0, "endpoint", &error)) {
        return make_error(env, error);
    }
    if (!find_boolean_in_prop_list(env, argv[4], "endpoint", &create_options.endpoint, true, &error)) {
        return make_error(env, error);
    }
    scalar_element_t start;
    scalar_element_t stop;
    scalar_element_t base;
    if (!parse_erl_scalar_term(env, argv[0], &start, &create_options.create_options, &error) ||
            !parse_erl_scalar_term(env, argv[1], &stop, &create_options.create_options, &error) ||
            !parse_erl_scalar_term(env, argv[3], &base, &create_options.create_options, &error)
            )
    {
        return make_error(env, error);
    }

    unsigned count = 0;
    if (!enif_get_uint(env, argv[2], &count)) {
        return make_error(env, ERR_ARG_BAD_COUNT);
    }
    nvector_t result;
    memset(&result, 0, sizeof(result));

    if (0==count) {
        result.array.dtype = DTAuto==create_options.create_options.dtype
                ? common_dtype(start.dtype, stop.dtype)
                : create_options.create_options.dtype;
        result.array.size0 = 1;
        result.array.size1 = 0;
        if (!enif_alloc_binary(0, &result.array.bin)) {
            return make_error(env, ERR_GEN_MEMORY_ALLOCATION);
        }
        return nvector_to_erl_record(env, &result, 0);
    }

    if (DTAuto==create_options.create_options.dtype) {
        double pstart = 0.0;
        double pstop = 0.0;
        double _Complex pcstart = 0.0;
        double _Complex pcstop = 0.0;
        _Bool start_single = false;
        _Bool stop_single = false;
        const _Bool complex_values =
                DTComplex==start.dtype || DTDoubleComplex==start.dtype ||
                DTComplex==stop.dtype || DTDoubleComplex==stop.dtype;
        if ( (DTSingle==start.dtype || DTDouble==start.dtype) && (DTSingle==base.dtype || DTDouble==base.dtype) ) {
            double dbase = double_value(&base, false);
            double dstart = double_value(&start, false);
            pstart = pow(dbase, dstart);
            start_single = maybe_single_precision(pstart);
        }
        else {
            double _Complex cbase = double_complex_value(&base);
            double _Complex cstart = double_complex_value(&start);
            pcstart = cpow(cbase, cstart);
            start_single = maybe_single_precision(creal(pcstart)) && maybe_single_precision(cimag(pcstart));
        }
        if ( (DTSingle==stop.dtype || DTDouble==stop.dtype) && (DTSingle==stop.dtype || DTDouble==stop.dtype) ) {
            double dbase = double_value(&base, false);
            double dstop = double_value(&stop, false);
            pstart = pow(dbase, dstop);
            stop_single = maybe_single_precision(pstop);
        }
        else {
            double _Complex cbase = double_complex_value(&base);
            double _Complex cstop = double_complex_value(&stop);
            pcstop = cpow(cbase, cstop);
            stop_single = maybe_single_precision(creal(pcstop)) && maybe_single_precision(cimag(pcstop));
        }
        if (complex_values) {
            result.array.dtype = start_single && stop_single
                    ? DTComplex : DTDoubleComplex;
        }
        else {
            result.array.dtype = start_single && stop_single
                    ? DTSingle : DTDouble;
        }
    }
    else {
        result.array.dtype = create_options.create_options.dtype;
    }

    result.array.size0 = 1;
    result.array.size1 = count;
    result.view.offset = 0;
    result.view.size = count;
    result.view.increment = 1;

    if (!narray_allocate_linspace(&result.array, &start, &stop, create_options.endpoint, &error)) {
        return make_error(env, error);
    }

    const double dbase = double_value(&base, false);
    const double _Complex zbase = double_complex_value(&base);
    const float fbase = (float) dbase;
    const float _Complex cbase = (float _Complex) zbase;

    // Convert linspace to logspace
    const size_t item_size = scalar_size(result.array.dtype);
    for (size_t i=0; i<count; ++i) {
        scalar_element_t src;
        src.dtype = result.array.dtype;
        memcpy(&src.value, result.array.bin.data + i*item_size, item_size);
        scalar_element_t dst;
        dst.dtype = src.dtype;
        switch (src.dtype) {
        case DTSingle:          dst.value.s = powf(fbase, src.value.s);     break;
        case DTDouble:          dst.value.d = pow(dbase, src.value.d);      break;
        case DTComplex:         dst.value.c = cpowf(cbase, src.value.c);    break;
        case DTDoubleComplex:   dst.value.z = cpow(zbase, src.value.z);     break;
        default:                memset(&dst.value, 0, sizeof(dst.value));   break;
        }
        memcpy(result.array.bin.data + i*item_size, &dst.value, item_size);
    }

    return nvector_to_erl_record(env, &result, 0);
}


static scalar_element_t
to_log10_value(const scalar_element_t v)
{
    scalar_element_t result;
    switch (v.dtype) {
    case DTSingle:
    case DTDouble:
    {
        double arg = DTSingle==v.dtype ? (double) v.value.s : v.value.d;
        double dlog = log10(arg);
        if (maybe_single_precision(dlog)) {
            result.dtype = DTSingle;
            result.value.s = (float) dlog;
        }
        else {
            result.dtype = DTDouble;
            result.value.d = dlog;
        }
        break;
    }

    case DTComplex:
    case DTDoubleComplex:
    {
        double _Complex arg = DTComplex==v.dtype ? (double _Complex) v.value.c : v.value.z;
        double _Complex zlog = clog10(arg);
        double zlog_re = creal(zlog);
        double zlog_im = cimag(zlog);
        result.dtype = maybe_single_precision(zlog_re) && maybe_single_precision(zlog_im)
                ? DTComplex : DTDoubleComplex;
        if (DTComplex==result.dtype) {
            result.value.c = (float _Complex) zlog;
        }
        else {
            result.value.z = zlog;
        }
        break;
    }
    default:
        memset(&result.value, 0, sizeof(result));
        result.dtype = DTSingle;
        break;
    }
    return result;
}


ERL_NIF_TERM
erl_nvector_geomspace(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (4!=argc) {
        return enif_make_badarg(env);
    }
    range_create_options_t create_options;

    const char* error = 0;
    if (!parse_create_options(env, argv[3], &create_options.create_options, 0, "endpoint", &error)) {
        return make_error(env, error);
    }
    if (!find_boolean_in_prop_list(env, argv[3], "endpoint", &create_options.endpoint, true, &error)) {
        return make_error(env, error);
    }
    scalar_element_t start;
    scalar_element_t stop;
    if (!parse_erl_scalar_term(env, argv[0], &start, &create_options.create_options, &error) ||
            !parse_erl_scalar_term(env, argv[1], &stop, &create_options.create_options, &error)
            )
    {
        return make_error(env, error);
    }

    const double start_re = double_value(&start, false);
    const double start_im = double_value(&start, true);
    if (0==start_re && 0==start_im) {
        return make_error(env, ERR_ARG_BAD_START);
    }

    const double stop_re = double_value(&stop, false);
    const double stop_im = double_value(&stop, true);
    if (0==stop_re && 0==stop_im) {
        return make_error(env, ERR_ARG_BAD_STOP);
    }

    uint8_t sign = 1;

    if (stop_re < 0.0 && start_re < 0.0) {
        sign = -1;
        scalar_element_t t = start;
        start = stop;
        stop = t;
    }

    start = to_log10_value(start);
    stop = to_log10_value(stop);

    unsigned count = 0;
    if (!enif_get_uint(env, argv[2], &count)) {
        return make_error(env, ERR_ARG_BAD_COUNT);
    }
    nvector_t result;
    memset(&result, 0, sizeof(result));

    if (0==count) {
        result.array.dtype = DTAuto==create_options.create_options.dtype
                ? common_dtype(start.dtype, stop.dtype)
                : create_options.create_options.dtype;
        result.array.size0 = 1;
        result.array.size1 = 0;
        result.view.offset = 0;
        result.view.size = 0;
        result.view.increment = 1;
        if (!enif_alloc_binary(0, &result.array.bin)) {
            return make_error(env, ERR_GEN_MEMORY_ALLOCATION);
        }
        return nvector_to_erl_record(env, &result, 0);
    }

    if (DTAuto==create_options.create_options.dtype) {
        double pstart = 0.0;
        double pstop = 0.0;
        double _Complex pcstart = 0.0;
        double _Complex pcstop = 0.0;
        _Bool start_single = false;
        _Bool stop_single = false;
        const _Bool complex_values =
                DTComplex==start.dtype || DTDoubleComplex==start.dtype ||
                DTComplex==stop.dtype || DTDoubleComplex==stop.dtype;
        if ( (DTSingle==start.dtype || DTDouble==start.dtype) ) {
            double dstart = double_value(&start, false);
            pstart = pow(10, dstart);
            start_single = maybe_single_precision(pstart);
        }
        else {
            double _Complex cstart = double_complex_value(&start);
            pcstart = cpow(10, cstart);
            start_single = maybe_single_precision(creal(pcstart)) && maybe_single_precision(cimag(pcstart));
        }
        if ( (DTSingle==stop.dtype || DTDouble==stop.dtype) && (DTSingle==stop.dtype || DTDouble==stop.dtype) ) {
            double dstop = double_value(&stop, false);
            pstart = pow(10, dstop);
            stop_single = maybe_single_precision(pstop);
        }
        else {
            double _Complex cstop = double_complex_value(&stop);
            pcstop = cpow(10, cstop);
            stop_single = maybe_single_precision(creal(pcstop)) && maybe_single_precision(cimag(pcstop));
        }
        if (complex_values) {
            result.array.dtype = start_single && stop_single
                    ? DTComplex : DTDoubleComplex;
        }
        else {
            result.array.dtype = start_single && stop_single
                    ? DTSingle : DTDouble;
        }
    }
    else {
        result.array.dtype = create_options.create_options.dtype;
    }

    result.array.size0 = 1;
    result.array.size1 = count;
    result.view.offset = 0;
    result.view.size = count;
    result.view.increment = 1;

    if (!narray_allocate_linspace(&result.array, &start, &stop, create_options.endpoint, &error)) {
        return make_error(env, error);
    }


    // Convert linspace to geomspace via logspace
    const size_t item_size = scalar_size(result.array.dtype);
    for (size_t i=0; i<count; ++i) {
        scalar_element_t src;
        src.dtype = result.array.dtype;
        memcpy(&src.value, result.array.bin.data + i*item_size, item_size);
        scalar_element_t dst;
        dst.dtype = src.dtype;
        switch (src.dtype) {
        case DTSingle:          dst.value.s = sign * powf(10, src.value.s);     break;
        case DTDouble:          dst.value.d = sign * pow(10, src.value.d);      break;
        case DTComplex:         dst.value.c = sign * cpowf(10, src.value.c);    break;
        case DTDoubleComplex:   dst.value.z = sign * cpow(10, src.value.z);     break;
        default:                memset(&dst.value, 0, sizeof(dst.value));   break;
        }
        memcpy(result.array.bin.data + i*item_size, &dst.value, item_size);
    }

    return nvector_to_erl_record(env, &result, 0);
}


ERL_NIF_TERM
erl_nvector_copy(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (2!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nvector_t vec;
    ERL_NIF_TERM bin_ref = 0;
    if (!parse_erl_vector(env, argv[0], &vec, &error, &bin_ref)) {
        return make_error(env, error);
    }
    create_options_t options;
    if (!parse_create_options(env, argv[1], &options, 0,  0, &error)) {
        return make_error(env, error);
    }
    dtype_t target_dtype = vec.array.dtype;
    if (DTAuto != options.dtype) {
        target_dtype = options.dtype;
    } else if (PSingle==options.precision && DTDouble==vec.array.dtype) {
        target_dtype = DTSingle;
    } else if (PSingle==options.precision && DTDoubleComplex==vec.array.dtype) {
        target_dtype = DTComplex;
    } else if (PDouble==options.precision && DTSingle==vec.array.dtype) {
        target_dtype = DTDouble;
    } else if (PDouble==options.precision && DTComplex==vec.array.dtype) {
        target_dtype = DTDoubleComplex;
    }
    nvector_t result;
    result.view.offset = 0;
    result.view.increment = 1;
    result.view.size = vec.view.size;
    result.array.dtype = target_dtype;
    result.array.size0 = 1;
    result.array.size1 = vec.view.size;

    const size_t item_size = scalar_size(target_dtype);
    if (target_dtype==vec.array.dtype && 1==vec.view.increment) {
        // Just return binary reference or make sub-binary, but not make full copy
        const size_t binary_items_count = vec.array.size0 * vec.array.size1;
        if (0==vec.view.offset && binary_items_count==vec.view.size) {
            // Unchanged binary
            return nvector_to_erl_record(env, &result, enif_make_copy(env, bin_ref));
        }
        else {
            // Continuous sub-binary
            const size_t bin_pos = item_size * vec.view.offset;
            const size_t bin_size = item_size * vec.view.size;
            return nvector_to_erl_record(env, &result, enif_make_sub_binary(env, bin_ref, bin_pos, bin_size));
        }
    }
    if (!enif_alloc_binary(vec.view.size*item_size, &result.array.bin)) {
        return make_error(env, ERR_GEN_MEMORY_ALLOCATION);
    }
    if (!narray_copy(&vec.array, result.array.bin.data, vec.view, target_dtype, &error)) {
        enif_release_binary(&result.array.bin);
        return make_error(env, error);
    }
    return nvector_to_erl_record(env, &result, 0);
}

ERL_NIF_TERM
erl_nvector_scale(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (3!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nvector_t vec;
    ERL_NIF_TERM bin_ref = 0;
    if (!parse_erl_vector(env, argv[0], &vec, &error, &bin_ref)) {
        return make_error(env, error);
    }
    create_options_t options;
    if (!parse_create_options(env, argv[2], &options, 0, 0, &error)) {
        return make_error(env, error);
    }
    scalar_element_t alpha;
    if (!parse_erl_scalar_term(env, argv[1], &alpha, &options, &error)) {
        return make_error(env, error);
    }
    dtype_t target_dtype = vec.array.dtype;
    if (DTAuto != options.dtype) {
        target_dtype = options.dtype;
    } else if (PSingle==options.precision && DTDouble==vec.array.dtype) {
        target_dtype = DTSingle;
    } else if (PSingle==options.precision && DTDoubleComplex==vec.array.dtype) {
        target_dtype = DTComplex;
    } else if (PDouble==options.precision && DTSingle==vec.array.dtype) {
        target_dtype = DTDouble;
    } else if (PDouble==options.precision && DTComplex==vec.array.dtype) {
        target_dtype = DTDoubleComplex;
    }
    nvector_t result;
    result.view.offset = 0;
    result.view.increment = 1;
    result.view.size = vec.view.size;
    result.array.dtype = target_dtype;
    result.array.size0 = 1;
    result.array.size1 = vec.view.size;

    const size_t item_size = scalar_size(target_dtype);
    if (!enif_alloc_binary(vec.view.size*item_size, &result.array.bin)) {
        return make_error(env, ERR_GEN_MEMORY_ALLOCATION);
    }
    if (!narray_scale(&vec.array, result.array.bin.data, vec.view, target_dtype, alpha, &error)) {
        enif_release_binary(&result.array.bin);
        return make_error(env, error);
    }
    return nvector_to_erl_record(env, &result, 0);
}




ERL_NIF_TERM
erl_nvector_axpy(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (4!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nvector_t x_vec, y_vec;
    if (!parse_erl_vector(env, argv[0], &y_vec, &error, NULL)) {
        return make_error(env, error);
    }
    if (!parse_erl_vector(env, argv[1], &x_vec, &error, NULL)) {
        return make_error(env, error);
    }
    create_options_t options;
    if (!parse_create_options(env, argv[3], &options, 0, 0, &error)) {
        return make_error(env, error);
    }
    scalar_element_t alpha;
    if (!parse_erl_scalar_term(env, argv[2], &alpha, &options, &error)) {
        return make_error(env, error);
    }
    dtype_t target_dtype = y_vec.array.dtype;
    if (DTAuto != options.dtype) {
        target_dtype = options.dtype;
    } else if (PSingle==options.precision && DTDouble==x_vec.array.dtype) {
        target_dtype = DTSingle;
    } else if (PSingle==options.precision && DTDoubleComplex==x_vec.array.dtype) {
        target_dtype = DTComplex;
    } else if (PDouble==options.precision && DTSingle==x_vec.array.dtype) {
        target_dtype = DTDouble;
    } else if (PDouble==options.precision && DTComplex==x_vec.array.dtype) {
        target_dtype = DTDoubleComplex;
    }
    nvector_t result;
    result.view.offset = 0;
    result.view.increment = 1;
    result.view.size = max_size(x_vec.view.size, y_vec.view.size);
    result.array.dtype = target_dtype;
    result.array.size0 = 1;
    result.array.size1 = result.view.size;

    const size_t item_size = scalar_size(target_dtype);
    if (!enif_alloc_binary(result.view.size*item_size, &result.array.bin)) {
        return make_error(env, ERR_GEN_MEMORY_ALLOCATION);
    }
    if (!narray_axpy(&x_vec.array,
                     &y_vec.array,
                     result.array.bin.data,
                     x_vec.view,
                     y_vec.view,
                     result.array.dtype,
                     alpha,
                     &error))
    {
        enif_release_binary(&result.array.bin);
        return make_error(env, error);
    }
    return nvector_to_erl_record(env, &result, 0);
}

ERL_NIF_TERM
erl_nvector_dot(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (3!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nvector_t x_vec, y_vec;
    if (!parse_erl_vector(env, argv[0], &y_vec, &error, NULL)) {
        return make_error(env, error);
    }
    if (!parse_erl_vector(env, argv[1], &x_vec, &error, NULL)) {
        return make_error(env, error);
    }
    _Bool conjuated = false;
    create_options_t options;
    if (!parse_create_options(env, argv[2], &options, 0, 0, &error)) {
        return make_error(env, error);
    }
    scalar_element_t res;
    if (DTAuto==options.dtype && PSingle!=options.precision)
        res.dtype = DTDouble;
    else
        res.dtype = DTSingle;
    if (!narray_dot(&x_vec.array,
                    &y_vec.array,
                    &res,
                    PSingle!=options.precision,
                    conjuated,
                    x_vec.view,
                    y_vec.view,
                    &error))
    {
        return make_error(env, error);
    }
    if (DTAuto!=options.dtype) {
        res.value = convert_type(res.value, res.dtype, options.dtype);
        res.dtype = options.dtype;
    }
    return make_scalar_value(env, res, CVTnoconvert);
}

ERL_NIF_TERM
erl_nvector_asum(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (2!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nvector_t vec;
    if (!parse_erl_vector(env, argv[0], &vec, &error, NULL)) {
        return make_error(env, error);
    }

    create_options_t options;
    if (!parse_create_options(env, argv[1], &options, 0, 0, &error)) {
        return make_error(env, error);
    }
    scalar_element_t res;
    if (!narray_asum(&vec.array, vec.view, &res, &error))
    {
        return make_error(env, error);
    }
    if (DTAuto!=options.dtype) {
        res.value = convert_type(res.value, res.dtype, options.dtype);
        res.dtype = options.dtype;
    }
    return make_scalar_value(env, res, CVTnoconvert);
}

ERL_NIF_TERM
erl_nvector_iamax_iamin(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (2!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nvector_t vec;
    if (!parse_erl_vector(env, argv[0], &vec, &error, NULL)) {
        return make_error(env, error);
    }
    char mode_atom[64]; memset(mode_atom, 0, sizeof(mode_atom));
    if (!enif_get_atom(env, argv[1], mode_atom, sizeof(mode_atom), ERL_NIF_LATIN1)) {
        return enif_make_badarg(env);
    }
    _Bool minMode;
    if (0==strncmp("iamin", mode_atom, sizeof(mode_atom))) {
        minMode = true;
    }
    else if (0==strncmp("iamax", mode_atom, sizeof(mode_atom))) {
        minMode = false;
    }
    else {
        return make_error(env, ERR_ARG_BAD_MIN_MAX_ATOM);
    }
    if (0==vec.view.size) {
        return make_atom(env, "undefined");
    }
    size_t result = 0;
    if (!narray_iamax_iamin(minMode, &vec.array, vec.view, &result, &error)) {
        return make_error(env, error);
    }
    return enif_make_uint64(env, result);
}
