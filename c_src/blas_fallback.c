#include "blas_fallback.h"
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static size_t
fallback_isamin(const int n, const float *x, const int incx)
{
    if (1==n || 0==n) { return 0; }
    float min_val = fabsf(x[0]);
    size_t min_index = 0;
    float cur_val;
    for (int i=1; i<n; ++i) {
        int real_index = i*incx;
        cur_val = fabsf(x[real_index]);
        if (cur_val < min_val) {
            min_index = i;
            min_val = cur_val;
        }
    }
    return min_index;
}

static size_t
fallback_idamin(const int n, const double *x, const int incx)
{
    if (1==n || 0==n) { return 0; }
    double min_val = fabs(x[0]);
    size_t min_index = 0;
    double cur_val;
    for (int i=1; i<n; ++i) {
        int real_index = i*incx;
        cur_val = fabs(x[real_index]);
        if (cur_val < min_val) {
            min_index = i;
            min_val = cur_val;
        }
    }
    return min_index;
}

static size_t
fallback_icamin(const int n, const float _Complex *x, const int incx)
{
    if (1==n || 0==n) { return 0; }
    float min_val = fabsf(crealf(x[0])) + fabsf(cimagf(x[0]));
    size_t min_index = 0;
    float cur_val;
    for (int i=1; i<n; ++i) {
        int real_index = i*incx;
        cur_val = fabsf(crealf(x[real_index])) + fabsf(cimagf(x[real_index]));
        if (cur_val < min_val) {
            min_index = i;
            min_val = cur_val;
        }
    }
    return min_index;
}

static size_t
fallback_izamin(const int n, const double _Complex *x, const int incx)
{
    if (1==n || 0==n) { return 0; }
    double min_val = fabs(creal(x[0])) + fabs(cimag(x[0]));
    size_t min_index = 0;
    double cur_val;
    for (int i=1; i<n; ++i) {
        int real_index = i*incx;
        cur_val = fabs(creal(x[real_index])) + fabs(cimag(x[real_index]));
        if (cur_val < min_val) {
            min_index = i;
            min_val = cur_val;
        }
    }
    return min_index;
}

typedef struct {
    const char * name;
    void * ptr;
} function_spec_t;

static const function_spec_t Functions_Available[] = {
    {"cblas_isamin", fallback_isamin},
    {"cblas_idamin", fallback_idamin},
    {"cblas_icamin", fallback_icamin},
    {"cblas_izamin", fallback_izamin},
    {NULL, NULL}
};


void*
resolve_fallback_function(const char *func_name)
{
    for (const function_spec_t *it=Functions_Available; it->name; ++it) {
        if (0==strcmp(func_name, it->name)) {
            return it->ptr;
        }
    }
    return 0;
}
