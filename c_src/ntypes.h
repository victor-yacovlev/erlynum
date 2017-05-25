#ifndef NTYPES_H
#define NTYPES_H

#include <complex.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

typedef enum precision {PAuto=0, PSingle=1, PDouble=2} precision_t;
typedef enum dtype {DTAuto=0, DTSingle=1, DTDouble=2, DTComplex=3, DTDoubleComplex=4} dtype_t;
typedef enum convert_option {CVTnoconvert=0, CVTinteger=1, CVTreal=2, CVTcomplex=3} convert_option_t;
typedef struct view_params {
    size_t          offset;
    size_t          size;
    int32_t         increment;
} view_params_t;
typedef unsigned char* raw_data_t;
typedef union {
    float           s;
    double          d;
    float _Complex  c;
    double _Complex z;
} scalar_value_t;

typedef struct scalar_element {
    dtype_t         dtype;
    scalar_value_t  value;
} scalar_element_t;

typedef struct create_options {
    precision_t     precision;
    dtype_t         dtype;
} create_options_t;

typedef struct range_create_options {
    create_options_t    create_options;
    _Bool               endpoint;
} range_create_options_t;

dtype_t
common_dtype(const dtype_t a, const dtype_t b);

dtype_t
common_dtype_for_elements(const size_t N, const scalar_element_t * items);

scalar_value_t
convert_type(const scalar_value_t src, const dtype_t from, const dtype_t to);

// return dtype is y_dtype
scalar_value_t
scalar_axpy(const scalar_value_t a, const dtype_t a_dtype,
            const scalar_value_t x, const dtype_t x_dtype,
            const scalar_value_t y, const dtype_t y_dtype);

void
convert_values(const size_t N, const void * src, void * dst,
               const dtype_t from, const dtype_t to);

double
double_value(const scalar_element_t *s, _Bool imaginary);

double _Complex
double_complex_value(const scalar_element_t *s);

size_t
scalar_size(const dtype_t dtype);

_Bool
maybe_single_precision(const double value);

#endif // NTYPES_H
