#include "ntypes.h"
#include <complex.h>
#include <float.h>
#include <string.h>

size_t
scalar_size(const dtype_t dtype)
{
    switch (dtype) {
    case DTSingle:          return sizeof(float);
    case DTDouble:          return sizeof(double);
    case DTComplex:         return sizeof(float _Complex);
    case DTDoubleComplex:   return sizeof(double _Complex);
    default:                return sizeof(void*);
    }
}

dtype_t
common_dtype(const dtype_t a, const dtype_t b)
{
    return a > b ? a : b;
}

dtype_t
common_dtype_for_elements(const size_t N, const scalar_element_t *items)
{
    dtype_t result = DTSingle;
    for (size_t i=0; i<N; ++i) {
        const scalar_element_t * current = items + i*sizeof(*items);
        const dtype_t current_dtype = current->dtype;
        if (DTAuto != current_dtype) {
            result = common_dtype(current_dtype, result);
        }
    }
    return result;
}

scalar_value_t convert_type(const scalar_value_t src, const dtype_t from, const dtype_t to)
{
    if (to==from)
        return src;
    scalar_value_t result;
    memset(&result, 0, sizeof(result));
    if (DTSingle==to) {
        switch (from) {
        case DTDouble:          result.s = (float) src.d;           break;
        case DTComplex:         result.s = crealf(src.c);           break;
        case DTDoubleComplex:   result.s = (float) creal(src.z);    break;
        default: break;
        }
    }
    else if (DTDouble==to) {
        switch (from) {
        case DTSingle:          result.d = src.s;                   break;
        case DTComplex:         result.d = crealf(src.c);           break;
        case DTDoubleComplex:   result.d = creal(src.z);            break;
        default: break;
        }
    }
    else if (DTComplex==to) {
        switch (from) {
        case DTSingle:          result.c = CMPLXF(src.s, 0.0);      break;
        case DTDouble:          result.c = CMPLXF((float)src.d, 0.0); break;
        case DTDoubleComplex:   result.c = CMPLXF(creal(src.z), cimag(src.z)); break;
        default: break;
        }
    }
    else if (DTDoubleComplex==to) {
        switch (from) {
        case DTSingle:          result.z = CMPLX(src.s, 0.0);       break;
        case DTDouble:          result.z = CMPLX(src.d, 0.0);       break;
        case DTComplex:         result.z = CMPLX(crealf(src.c), cimagf(src.c)); break;
        default: break;
        }
    }
    return result;
}

void
convert_values(const size_t N, const void * src, void * dst,
               const dtype_t from, const dtype_t to)
{
    if (from==to) {
        memcpy(dst, src, N * scalar_size(from));
        return;
    }
    scalar_value_t x, y;
    const size_t x_size = scalar_size(from);
    const size_t y_size = scalar_size(to);
    for (size_t i=0; i<N; ++i) {
        memcpy(&x, src + i*x_size, x_size);
        y = convert_type(x, from, to);
        memcpy(dst + i*y_size, &y, y_size);
    }
}

double
double_value(const scalar_element_t *s, _Bool imaginary)
{
    switch (s->dtype) {
    case DTSingle:          return imaginary ? 0.0 : s->value.s;
    case DTDouble:          return imaginary ? 0.0 : s->value.d;
    case DTComplex:         return imaginary ? cimagf(s->value.c) : creal(s->value.c);
    case DTDoubleComplex:   return imaginary ? cimag(s->value.z) : creal(s->value.z);
    default:                return 0.0;
    }
}

double _Complex
double_complex_value(const scalar_element_t *s)
{
    switch (s->dtype) {
    case DTSingle:          return CMPLX(s->value.s, 0);
    case DTDouble:          return CMPLX(s->value.d, 0);
    case DTComplex:         return CMPLX(crealf(s->value.c), cimagf(s->value.c));
    case DTDoubleComplex:   return s->value.z;
    default:                return CMPLX(0, 0);
    }
}

_Bool
maybe_single_precision(const double value)
{
    // 0 is 0
    if (0==value)
        return true;

    // Too big absolute value for float
    if (value < FLT_MIN || value > FLT_MAX)
        return false;

    // Check for precision loss
    volatile const float f = (float) value;
    volatile const double d = (double) f;
    volatile const _Bool matches = d==value;
    return matches;
}

scalar_value_t
scalar_axpy(const scalar_value_t a,
            const dtype_t a_dtype,
            const scalar_value_t x,
            const dtype_t x_dtype,
            const scalar_value_t y,
            const dtype_t y_dtype)
{
    const scalar_value_t norm_a = convert_type(a, a_dtype, y_dtype);
    scalar_value_t res = convert_type(x, x_dtype, y_dtype);
    if (DTSingle==y_dtype) {
        res.s *= norm_a.s;
        res.s += y.s;
    }
    else if (DTDouble==y_dtype) {
        res.d *= norm_a.d;
        res.d += y.d;
    }
    else if (DTComplex==y_dtype) {
        res.c *= norm_a.c;
        res.c += y.c;
    }
    else if (DTDoubleComplex==y_dtype) {
        res.z *= norm_a.z;
        res.z += y.z;
    }
    return res;
}
