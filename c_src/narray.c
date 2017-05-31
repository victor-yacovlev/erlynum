#include "narray.h"

#include "error_atoms.h"
#include "nblas.h"

#include <erl_nif.h>

#include <string.h>

_Bool
narray_allocate_range(narray_t *array,
                      const scalar_element_t *start,
                      const scalar_element_t *stop,
                      const scalar_element_t *step,
                      const _Bool endpoint,
                      const char* *error)
{
    double start_re = double_value(start, 0);
    double stop_re = double_value(stop, 0);
    double step_re = double_value(step, 0);

    double start_im = double_value(start, 1);
    double stop_im = double_value(stop, 1);
    double step_im = double_value(step, 1);

    if (0.0 == step_re && 0.0 == step_im) {
        *error = ERR_ARG_BAD_STEP;
        return false;
    }
    const double sign_re = step_re < 0.0 ? -1.0 : 1.0;
    const double sign_im = step_im < 0.0 ? -1.0 : 1.0;

    double range_re = stop_re - start_re;
    if (range_re < 0.0) {
        *error = ERR_ARG_BAD_RANGE;
        return false;
    }
    double range_im = stop_im - start_im;
    if (range_im < 0.0) {
        *error = ERR_ARG_BAD_RANGE;
        return false;
    }
    array->size0 = 1;
    array->size1 = 1;
    size_t steps_re = 0;
    size_t steps_im = 0;
    if (0.0!=step_re) {
        steps_re = (size_t)(sign_re * range_re / step_re);
    }
    if (0.0 != step_im) {
        steps_im = (size_t)(sign_im * range_im / step_im);
    }
    array->size1 = steps_re > steps_im ? steps_re : steps_im;

    const size_t alloc_size = endpoint ? (array->size1+1) : (array->size1);

    const size_t item_size = scalar_size(array->dtype);
    if (!enif_alloc_binary(item_size*alloc_size, &array->bin)) {
        *error = ERR_GEN_MEMORY_ALLOCATION;
        return false;
    }
    double current_re;
    double current_im;
    scalar_value_t current_value;
    for (size_t i=0; i<array->size1; ++i) {
        current_re = start_re + i*sign_re*step_re;
        current_im = start_im + i*sign_im*step_im;
        switch (array->dtype) {
        case DTSingle:          current_value.s = (float) current_re; break;
        case DTDouble:          current_value.d = current_re; break;
        case DTComplex:         current_value.c = CMPLXF((float)current_re, (float)current_im); break;
        case DTDoubleComplex:   current_value.z = CMPLX(current_re, current_im); break;
        default:                memset(&current_value, 0, sizeof(current_value)); break;
        }
        void * dst = array->bin.data + i*item_size;
        memcpy(dst, &current_value, item_size);
    }
    if (endpoint) {
        current_re = sign_re * stop_re;
        current_im = sign_im * stop_im;
        switch (array->dtype) {
        case DTSingle:          current_value.s = (float) current_re; break;
        case DTDouble:          current_value.d = current_re; break;
        case DTComplex:         current_value.c = CMPLXF((float)current_re, (float)current_im); break;
        case DTDoubleComplex:   current_value.z = CMPLX(current_re, current_im); break;
        default:                memset(&current_value, 0, sizeof(current_value)); break;
        }
        void * dst = array->bin.data + array->size1*item_size;
        memcpy(dst, &current_value, item_size);
    }
    array->size1 = alloc_size;
    return true;
}

_Bool
narray_allocate_linspace(narray_t * arr,
                       const scalar_element_t *start,
                       const scalar_element_t *stop,
                       const _Bool endpoint,
                       const char* *error)
{
    const size_t item_size = scalar_size(arr->dtype);
    const size_t N = arr->size0 * arr->size1;
    if (!enif_alloc_binary(item_size*N, &arr->bin)) {
        *error = ERR_GEN_MEMORY_ALLOCATION;
        return false;
    }
    const double start_re = double_value(start, false);
    const double start_im = double_value(start, true);
    const double stop_re = double_value(stop, false);
    const double stop_im = double_value(stop, true);
    const double range_re = stop_re - start_re;
    if (range_re < 0.0) {
        *error = ERR_ARG_BAD_RANGE;
        return false;
    }
    const double range_im = stop_im - start_im;
    if (range_im < 0.0) {
        *error = ERR_ARG_BAD_RANGE;
        return false;
    }
    const double step_re = endpoint ? range_re / (N-1) : range_re / N;
    const double step_im = endpoint ? range_im / (N-1) : range_im / N;
    const double sign_re = step_re < 0.0 ? -1.0 : 1.0;
    const double sign_im = step_im < 0.0 ? -1.0 : 1.0;

    double current_re;
    double current_im;
    scalar_value_t current_value;

    for (size_t i=0; i<N; ++i) {
        current_re = start_re + i*sign_re*step_re;
        current_im = start_im + i*sign_im*step_im;
        switch (arr->dtype) {
        case DTSingle:          current_value.s = (float) current_re; break;
        case DTDouble:          current_value.d = current_re; break;
        case DTComplex:         current_value.c = CMPLXF((float)current_re, (float)current_im); break;
        case DTDoubleComplex:   current_value.z = CMPLX(current_re, current_im); break;
        default:                memset(&current_value, 0, sizeof(current_value)); break;
        }
        void * dst = arr->bin.data + i*item_size;
        memcpy(dst, &current_value, item_size);
    }
    return true;
}

_Bool
narray_copy(const narray_t *x, void *y,
            const view_params_t x_range,
            const dtype_t y_dtype,
            const char **error)
{
    const size_t x_size = scalar_size(x->dtype);
    if (x->dtype==y_dtype) {
        // No data conversion, so can use BLAS copy routine
        return nblas_copy(x_range.size,
                          x->bin.data + x_size*x_range.offset,
                          x_range.increment,
                          y, 1,
                          y_dtype, error);
    }
    else {
        scalar_value_t x_val, y_val;
        const size_t y_size = scalar_size(y_dtype);
        for (size_t i=0; i<x_range.size; ++i) {
            const size_t real_index = x_range.offset + i * x_range.increment;
            memcpy(&x_val, x->bin.data + real_index*x_size, x_size);
            y_val = convert_type(x_val, x->dtype, y_dtype);
            memcpy(y + i*y_size, &y_val, y_size);
        }
        return true;
    }
}

_Bool
narray_scale(const narray_t *x, void *y,
             const view_params_t x_range,
             const dtype_t y_dtype,
             const scalar_element_t a,
             const char **error)
{
    const scalar_value_t norm_alpha = convert_type(a.value, a.dtype, y_dtype);
    if (!narray_copy(x, y, x_range, y_dtype, error))
        return false;
    return nblas_scal(x_range.size, norm_alpha, y, 1, y_dtype, error);
}



_Bool
narray_axpy(const narray_t *x,
            const narray_t *y,
            void *res,
            const view_params_t x_range,
            const view_params_t y_range,
            const dtype_t r_dtype,
            const scalar_element_t a,
            const char **error)
{
    const narray_t * real_x = 0;
    const narray_t * real_y = 0;
    view_params_t real_x_range, real_y_range;
    _Bool scale_result = false;
    view_params_t scale_range;
    // Y should be a vector of greater size
    if (y_range.size > x_range.size) {
        real_x = x;
        real_y = y;
        real_x_range = x_range;
        real_y_range = y_range;
        scale_result = false;
    }
    else {
        real_x = y;
        real_y = x;
        real_x_range = y_range;
        real_y_range = x_range;
        // Calculate a*x+0 for the rest:
        // just scale in case if x swapped with y
        scale_result = true;
        scale_range.size = x_range.size - y_range.size;
        scale_range.increment = x_range.increment;
        scale_range.offset = x_range.increment * y_range.size;
    }
    if (!narray_copy(real_y, res, real_y_range, r_dtype, error)) {
        return false;
    }
    const size_t x_size = scalar_size(x->dtype);
    const scalar_value_t norm_alpha = convert_type(a.value, a.dtype, r_dtype);
    if (x->dtype==r_dtype) {
        if (!nblas_axpy(real_x_range.size,
                          norm_alpha,
                          real_x->bin.data + x_size*real_x_range.offset,
                          real_x_range.increment,
                          res,
                          1,
                          r_dtype,
                          error))
        {
            return false;
        }
    }
    else {
        scalar_value_t x_val, r_val;
        const size_t r_size = scalar_size(r_dtype);
        for (size_t i=0; i<real_x_range.size; ++i) {
            const size_t real_index = real_x_range.offset + i * real_x_range.increment;
            memcpy(&x_val, real_x->bin.data + real_index*x_size, x_size);
            memcpy(&r_val, res + i*r_size, r_size);
            r_val = scalar_axpy(norm_alpha, r_dtype, x_val, x->dtype, r_val, r_dtype);
            memcpy(res + i*r_size, &r_val, r_size);
        }
    }
    if (scale_result) {
        const size_t r_size = scalar_size(r_dtype);
        return nblas_scal(scale_range.size,
                          norm_alpha,
                          res + scale_range.offset*r_size,
                          scale_range.increment,
                          r_dtype,
                          error);
    }
    else {
        return true;
    }
}

_Bool
narray_dot(const narray_t *x, const narray_t *y,
           scalar_element_t *res,
           _Bool double_precision,
           _Bool conjugated_complex,
           const view_params_t x_range,
           const view_params_t y_range,
           const char **error)
{
    const size_t x_n = x_range.size;
    const size_t y_n = y_range.size;
    const size_t c_n = x_n < y_n ? x_n : y_n;
    void *real_x = 0, *real_y = 0;
    const dtype_t c_dtype = common_dtype(x->dtype, y->dtype);
    const size_t item_size = scalar_size(c_dtype);
    view_params_t real_x_range = x_range;
    view_params_t real_y_range = y_range;
    real_x_range.size = real_y_range.size = c_n;
    if (c_dtype==x->dtype) {
        real_x = x->bin.data;
    }
    else {
        real_x = enif_alloc(item_size*c_n);
        if (!real_x) {
            *error = ERR_GEN_MEMORY_ALLOCATION;
            return false;
        }
        if (!narray_copy(x, real_x, real_x_range, c_dtype, error)) {
            enif_free(real_x);
            return false;
        }
        real_x_range.increment = 1;
        real_x_range.offset = 0;
    }
    if (c_dtype==y->dtype) {
        real_y = y->bin.data;
    }
    else {
        real_y = enif_alloc(item_size*c_n);
        if (!real_y) {
            *error = ERR_GEN_MEMORY_ALLOCATION;
            if (real_x != x->bin.data) {
                enif_free(real_x);
            }
            return false;
        }
        if (!narray_copy(y, real_y, real_y_range, c_dtype, error)) {
            enif_free(real_y);
            if (real_x != x->bin.data) {
                enif_free(real_x);
            }
            return false;
        }
        real_y_range.increment = 1;
        real_y_range.offset = 0;
    }

    _Bool status = false;
    scalar_element_t t;
    t.dtype = c_dtype;
    if (double_precision && DTSingle==c_dtype) {
        status = nblas_sdot(c_n,
                            (float*)(real_x + real_x_range.offset*sizeof(float)),
                            real_x_range.increment,
                            (float*)(real_y + real_y_range.offset*sizeof(float)),
                            real_y_range.increment,
                            &res->value, res->dtype, error);
    }
    else if (DTComplex==c_dtype || DTDoubleComplex==c_dtype) {
        status = nblas_dotu_dotc(c_n,
                                 real_x + real_x_range.offset*item_size,
                                 real_x_range.increment,
                                 real_y + real_y_range.offset*item_size,
                                 real_y_range.increment,
                                 conjugated_complex, &t.value, t.dtype, error);
    }
    else if (DTComplex!=c_dtype && DTDoubleComplex!=c_dtype) {
        status = nblas_dot(c_n,
                           real_x + real_x_range.offset*item_size,
                           real_x_range.increment,
                           real_y + real_y_range.offset*item_size,
                           real_y_range.increment,
                           &t.value, t.dtype, error);

    }
    if (real_x != x->bin.data) {
        enif_free(real_x);
    }
    if (real_y != x->bin.data) {
        enif_free(real_y);
    }
    if (!status) {
        res->value = convert_type(t.value, t.dtype, res->dtype);
    }
    return status;
}

_Bool
narray_asum(const narray_t *x,
            const view_params_t x_range,
            scalar_element_t *res,
            const char **error)
{
    if (DTSingle==x->dtype || DTComplex==x->dtype)
        res->dtype = DTSingle;
    else
        res->dtype = DTDouble;
    const size_t item_size = scalar_size(x->dtype);
    const void * data_ptr = x->bin.data + item_size * x_range.offset;
    return nblas_asum(x_range.size, data_ptr, x_range.increment, &res->value, x->dtype, error);
}
