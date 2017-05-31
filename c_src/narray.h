#ifndef NARRAY_H
#define NARRAY_H

#include "ntypes.h"
#include <stdbool.h>
#include <stdint.h>

#include <erl_nif.h>


#ifdef __cplusplus
extern "C" {
#endif

typedef struct array {
    dtype_t         dtype;
    size_t          size0;
    size_t          size1;
    ErlNifBinary    bin;
} narray_t;


_Bool
narray_allocate_range(narray_t * arrray,
                      const scalar_element_t *start,
                      const scalar_element_t *stop,
                      const scalar_element_t *step,
                      const _Bool endpoint,
                      const char* *error);

_Bool
narray_allocate_linspace(narray_t *arr,
                       const scalar_element_t *start,
                       const scalar_element_t *stop,
                       const _Bool endpoint,
                       const char* *error);

_Bool
narray_copy(const narray_t * x,
            void * y,
            const view_params_t x_range,
            const dtype_t y_dtype,
            const char* *error);

_Bool
narray_scale(const narray_t * x,
             void * y,
             const view_params_t x_range,
             const dtype_t y_dtype,
             const scalar_element_t a,
             const char* *error);

_Bool
narray_axpy(const narray_t *x,
            const narray_t *y,
            void * res,
            const view_params_t x_range,
            const view_params_t y_range,
            const dtype_t r_dtype,
            const scalar_element_t a,
            const char* *error);

_Bool
narray_dot(const narray_t *x,
           const narray_t *y,
           scalar_element_t * res,
           _Bool double_precision,
           _Bool conjugated_complex,
           const view_params_t x_range,
           const view_params_t y_range,
           const char* *error);

_Bool
narray_asum(const narray_t *x,
            const view_params_t x_range,
            scalar_element_t * res,
            const char* *error);

_Bool
narray_iamax_iamin(_Bool minMode,
                   const narray_t *x,
                   const view_params_t x_range,
                   size_t *res,
                   const char* *error);

_Bool
narray_nrm2(const narray_t *x,
            const view_params_t x_range,
            scalar_element_t * res,
            const char* *error);


#ifdef __cplusplus
} // extern "C"
#endif

#endif // NARRAY_H
