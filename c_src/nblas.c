#include "nblas.h"

#include "backend_loader.h"
#include "error_atoms.h"
#include "erl_util.h"
#include "nvector.h"
#include "ntypes.h"

#include <erl_nif.h>
#include <string.h>


_Bool
nblas_check_range(const size_t size,
                  const size_t n_iter,
                  const size_t incr)
{
    const size_t required = 1 + (n_iter-1) * incr;
    return size >= required;
}


_Bool
nblas_copy(const int n, const void *x, const int incx,
           void *y, const int incy,
           const dtype_t dtype, const char* *error)
{
    char blas_func_name[128];
    nblas_constuct_func_name("copy", 1, &dtype, sizeof(blas_func_name), blas_func_name);
    static void *func = 0;
    if (!func)
        func = resolve_blas_function(blas_func_name);
    if (!func) {
        *error = ERR_BLAS_NOTFOUND_COPY;
        return false;
    }
    void (*f_ptr)(const int, const void*, const int, void*, const int);
    f_ptr = func;
    (*f_ptr)(n, x, incx, y, incy);
    return true;
}

_Bool
nblas_scal(const int n, const scalar_value_t a,
           void *x, const int incx,
           const dtype_t dtype, const char **error)
{
    char blas_func_name[128];
    nblas_constuct_func_name("scal", 1, &dtype, sizeof(blas_func_name), blas_func_name);
    static void *func = 0;
    if (!func)
        func = resolve_blas_function(blas_func_name);
    if (!func) {
        *error = ERR_BLAS_NOTFOUND_COPY;
        return false;
    }
    if (DTSingle==dtype) {
        const float fa = a.s;
        void (*f_ptr)(const int, const float, float*, const int);
        f_ptr = func;
        (*f_ptr)(n, fa, x, incx);
    }
    else if (DTDouble==dtype) {
        const double da = a.d;
        void (*f_ptr)(const int, const double, double*, const int);
        f_ptr = func;
        (*f_ptr)(n, da, x, incx);
    }
    else if (DTComplex==dtype) {
        const float _Complex ca = a.c;
        void (*f_ptr)(const int, const void*, void*, const int);
        f_ptr = func;
        (*f_ptr)(n, &ca, x, incx);
    }
    else if (DTDoubleComplex==dtype) {
        const double _Complex za = a.z;
        void (*f_ptr)(const int, const void*, void*, const int);
        f_ptr = func;
        (*f_ptr)(n, &za, x, incx);
    }
    return true;
}

_Bool
nblas_axpy(const int n, const scalar_value_t a,
           const void *x, const int incx,
           void *y, const int incy,
           const dtype_t dtype, const char **error)
{
    char blas_func_name[128];
    nblas_constuct_func_name("axpy", 1, &dtype, sizeof(blas_func_name), blas_func_name);
    static void *func = 0;
    if (!func)
        func = resolve_blas_function(blas_func_name);
    if (!func) {
        *error = ERR_BLAS_NOTFOUND_AXPY;
        return false;
    }
    if (DTSingle==dtype) {
        const float fa = a.s;
        void (*f_ptr)(const int, const float, const float*, const int, float*, const int);
        f_ptr = func;
        (*f_ptr)(n, fa, x, incx, y, incy);
    }
    else if (DTDouble==dtype) {
        const double da = a.d;
        void (*f_ptr)(const int, const double, const double*, const int, double*, const int);
        f_ptr = func;
        (*f_ptr)(n, da, x, incx, y, incy);
    }
    else if (DTComplex==dtype) {
        const float _Complex ca = a.c;
        void (*f_ptr)(const int, const void*, const void*, const int, void*, const int);
        f_ptr = func;
        (*f_ptr)(n, &ca, x, incx, y, incy);
    }
    else if (DTDoubleComplex==dtype) {
        const double _Complex za = a.z;
        void (*f_ptr)(const int, const void*, const void*, const int, void*, const int);
        f_ptr = func;
        (*f_ptr)(n, &za, x, incx, y, incy);
    }
    return true;
}

_Bool
nblas_dot(const int n,
          const void *x, const int incx,
          const void *y, const int incy,
          void * out,
          const dtype_t dtype, const char **error)
{
    char blas_func_name[128];
    nblas_constuct_func_name("dot", 1, &dtype, sizeof(blas_func_name), blas_func_name);
    static void *func = 0;
    if (!func)
        func = resolve_blas_function(blas_func_name);
    if (!func) {
        *error = ERR_BLAS_NOTFOUND_DOT;
        return false;
    }
    if (DTSingle==dtype) {
        float (*f_ptr)(const int, const float*, const int, const float*, const int);
        f_ptr = func;
        float res = (*f_ptr)(n, x, incx, y, incy);
        memcpy(out, &res, sizeof(res));
    }
    else if (DTDouble==dtype) {
        double (*f_ptr)(const int, const double*, const int, const double*, const int);
        f_ptr = func;
        double res = (*f_ptr)(n, x, incx, y, incy);
        memcpy(out, &res, sizeof(res));
    }
    else {
        *error = ERR_BLAS_NOTSUPPORT_DTYPE;
        return false;
    }
    return true;
}



_Bool
nblas_sdot(const int n,
           const float *sx, const int incx,
           const float *sy, const int incy,
           void *out,
           const dtype_t dtype, const char **error)
{
    char blas_func_name[128];
    if (DTSingle==dtype) {
        nblas_constuct_func_name("dsdot", 1, &dtype, sizeof(blas_func_name), blas_func_name);
    }
    else if (DTDouble==dtype) {
        nblas_constuct_func_name("sdot", 1, &dtype, sizeof(blas_func_name), blas_func_name);
    }
    else {
        *error = ERR_BLAS_NOTSUPPORT_DTYPE;
        return false;
    }

    static void *func = 0;
    if (!func)
        func = resolve_blas_function(blas_func_name);
    if (!func) {
        *error = ERR_BLAS_NOTFOUND_SDOT;
        return false;
    }
    if (DTSingle==dtype) {
        float (*f_ptr)(const int, const float, const float*, const int, const float*, const int);
        f_ptr = func;
        float res = (*f_ptr)(n, 0, sx, incx, sy, incy);
        memcpy(out, &res, sizeof(res));
    }
    else if (DTDouble==dtype) {
        double (*f_ptr)(const int, const float*, const int, const float*, const int);
        f_ptr = func;
        double res = (*f_ptr)(n, sx, incx, sy, incy);
        memcpy(out, &res, sizeof(res));
    }
    else {
        *error = ERR_BLAS_NOTSUPPORT_DTYPE;
        return false;
    }
    return true;
}


void
nblas_constuct_func_name(const char *main_name,
                         const size_t dtypes_count,
                         const dtype_t *dtypes,
                         const size_t out_size,
                         char *out)
{
    static const char* CBlasPrefix = "cblas_";
    static const size_t CBlasPrefixLen = 6;
    memset(out, 0, out_size);
    for (size_t i=0; i<dtypes_count; ++i) {
        char p;
        const dtype_t dt = dtypes[i];
        switch (dt) {
        case DTSingle:              p = 's'; break;
        case DTDouble:              p = 'd'; break;
        case DTComplex:             p = 'c'; break;
        case DTDoubleComplex:       p = 'z'; break;
        default:                    p = '_'; break;
        }
        out[CBlasPrefixLen+i] = p;
    }
    strncpy(out, CBlasPrefix, CBlasPrefixLen);
    strncpy(out+CBlasPrefixLen+dtypes_count, main_name, out_size-dtypes_count-CBlasPrefixLen);
}



