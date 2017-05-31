#ifndef NBLAS_H
#define NBLAS_H

#include "nvector.h"
#include "ntypes.h"

#include <erl_nif.h>
#include <stdbool.h>

void
nblas_constuct_func_name(
        const char * main_name,
        const size_t dtypes_count,
        const dtype_t * dtypes,
        const size_t out_size,
        char * out
        );

_Bool
nblas_check_range(const size_t size,
                  const size_t n_iter,
                  const size_t incr);

_Bool
nblas_copy(const int n, const void *x, const int incx,
           void *y, const int incy,
           const dtype_t dtype, const char* *error);

_Bool
nblas_scal(const int n, const scalar_value_t a,
           void *x, const int incx,
           const dtype_t dtype, const char* *error);

_Bool
nblas_axpy(const int n, const scalar_value_t a,
           const void *x, const int incx,
           void *y, const int incy,
           const dtype_t dtype, const char* *error);

_Bool
nblas_dot(const int n,
          const void *x, const int incx,
          const void *y, const int incy,
          void * out,
          const dtype_t dtype, const char* *error);

_Bool
nblas_sdot(const int n,
           const float *sx, const int incx,
           const float *sy, const int incy,
           void * out,
           const dtype_t dtype, const char* *error);

_Bool
nblas_dotu_dotc(const int n,
                const void *x, const int incx,
                const void *y, const int incy,
                _Bool conjugated,
                void * out,
                const dtype_t dtype, const char* *error);

_Bool
nblas_asum(const int n,
           const void *x, const int incx,
           void * out,
           const dtype_t dtype, const char* *error);

_Bool
nblas_iamax_iamin(const int n,
                  _Bool minMode,
                  const void *x, const int incx,
                  size_t * out,
                  const dtype_t dtype, const char* *error);

#endif // NBLAS_H
