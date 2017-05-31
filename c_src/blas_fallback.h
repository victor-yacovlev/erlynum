#ifndef BLAS_FALLBACK_H
#define BLAS_FALLBACK_H

#ifdef __cplusplus
extern "C" {
#endif

void*
resolve_fallback_function(const char * func_name);

#ifdef __cplusplus
}
#endif

#endif // BLAS_FALLBACK_H
