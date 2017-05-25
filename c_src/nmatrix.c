#include "erl_util.h"
#include "error_atoms.h"
#include "ntypes.h"
#include "nmatrix.h"

#include <erl_nif.h>
#include <string.h>

ERL_NIF_TERM
erl_nmatrix_full(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (3!=argc) {
        return enif_make_badarg(env);
    }

    const char * error = 0;
    create_options_t options;
    if (!parse_create_options(env, argv[2], &options, 0, &error)) {
        return make_error(env, error);
    }

    nmatrix_t result;   memset(&result, 0, sizeof(result));
    scalar_element_t fill_value;

    if (!parse_erl_scalar_term(env, argv[1], &fill_value, &options, &error)) {
        return make_error(env, error);
    }
    result.array.dtype = fill_value.dtype;
    size_t shape[2];

    if (!parse_erl_shape(env, argv[0], 2, shape, &error)) {
        return make_error(env, error);
    }

    result.array.size0 = shape[0];
    result.array.size1 = shape[1];
    result.view0.offset = 0;
    result.view0.size = shape[0];
    result.view0.increment = (int32_t) shape[1];
    result.view1.offset = 0;
    result.view1.size = shape[1];
    result.view1.increment = 1;

    const size_t item_size = scalar_size(result.array.dtype);
    const size_t bin_size = item_size * shape[0] * shape[1];
    if (!enif_alloc_binary(bin_size, &result.array.bin)) {
        return make_error(env, ERR_GEN_MEMORY_ALLOCATION);
    }
    void * p = result.array.bin.data;
    for (size_t i=0; i<shape[0]*shape[1]; ++i) {
        memcpy(p, &fill_value.value, item_size);
        p += item_size;
    }

    return nmatrix_to_erl_record(env, &result, 0);
}

ERL_NIF_TERM
erl_nmatrix_eye(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (4!=argc) {
        return enif_make_badarg(env);
    }
    const char * error = 0;
    create_options_t options;
    if (!parse_create_options(env, argv[2], &options, 0, &error)) {
        return make_error(env, error);
    }

    nmatrix_t result;   memset(&result, 0, sizeof(result));
    scalar_element_t fill_value;
    if (!parse_erl_scalar_term(env, argv[1], &fill_value, &options, &error)) {
        return make_error(env, error);
    }
    result.array.dtype = fill_value.dtype;

    size_t shape[2];

    if (!parse_erl_shape(env, argv[0], 2, shape, &error)) {
        return make_error(env, error);
    }
    result.array.size0 = shape[0];
    result.array.size1 = shape[1];
    result.view0.offset = 0;
    result.view0.size = shape[0];
    result.view0.increment = (int32_t) shape[1];
    result.view1.offset = 0;
    result.view1.size = shape[1];
    result.view1.increment = 1;
    int diag_number = 0;
    if (!enif_get_int(env, argv[2], &diag_number)) {
        return make_error(env, ERR_ARG_BAD_DIAGONAL);
    }
    const size_t item_size = scalar_size(result.array.dtype);
    const size_t bin_size = item_size * shape[0] * shape[1];
    if (!enif_alloc_binary(bin_size, &result.array.bin)) {
        return make_error(env, ERR_GEN_MEMORY_ALLOCATION);
    }
    void * p = result.array.bin.data;
    memset(p, 0, bin_size);
    size_t x, y, plain_index;
    const size_t rows = shape[0];
    const size_t cols = shape[1];
    for (x=diag_number, y=0; x<cols && y<rows; ++x, ++y) {
        plain_index = y*cols + x;
        memcpy(p + plain_index*item_size, &fill_value.value, item_size);
    }
    return nmatrix_to_erl_record(env, &result, 0);
}

ERL_NIF_TERM
nmatrix_to_erl_record(ErlNifEnv *env, nmatrix_t *mat, const ERL_NIF_TERM bin_ref)
{
    const ERL_NIF_TERM header   = make_atom(env, "nmatrix");
    const ERL_NIF_TERM view0    = make_view_params(env, &mat->view0);
    const ERL_NIF_TERM view1    = make_view_params(env, &mat->view1);
    const ERL_NIF_TERM shape    = make_shape(env, mat->array.size0, mat->array.size1);
    const ERL_NIF_TERM dtype    = make_dtype(env, mat->array.dtype);
    const ERL_NIF_TERM data     = bin_ref ? bin_ref : enif_make_binary(env, &mat->array.bin);
    return enif_make_tuple5(env,
                            header,
                            enif_make_tuple2(env, view0, view1),
                            shape,
                            dtype,
                            data);
}

ERL_NIF_TERM
erl_nmatrix_to_list(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (2!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nmatrix_t mat;
    if (!parse_erl_matrix(env, argv[0], &mat, &error, NULL)) {
        return make_error(env, error);
    }
    convert_option_t convert;
    if (!parse_convert_option(env, argv[1], &convert, &error)) {
        return make_error(env, error);
    }
    const size_t rows_count = mat.view0.size;
    const size_t cols_count = mat.view1.size;
    const size_t item_size = scalar_size(mat.array.dtype);
    ERL_NIF_TERM * list_items = enif_alloc(rows_count * sizeof(ERL_NIF_TERM));
    for (size_t row_no=0; row_no<rows_count; ++row_no) {
        ERL_NIF_TERM * row_list_items = enif_alloc(cols_count * sizeof(ERL_NIF_TERM));
        for (size_t col_no=0; col_no<cols_count; ++col_no) {
            const int32_t view0_offset_index =
                    mat.view0.offset + row_no * mat.view0.increment;
            const int32_t view1_offset_index =
                    mat.view1.offset + col_no * mat.view1.increment;
            const size_t plain_index = view0_offset_index + view1_offset_index;
            void * src = mat.array.bin.data + plain_index * item_size;
            scalar_element_t element;
            element.dtype = mat.array.dtype;
            memcpy(&element.value, src, item_size);
            row_list_items[col_no] = make_scalar_value(env, element, convert);
        }
        list_items[row_no] = enif_make_list_from_array(env, row_list_items, cols_count);
        enif_free(row_list_items);
    }
    ERL_NIF_TERM result = enif_make_list_from_array(env, list_items, rows_count);
    enif_free(list_items);
    return result;
}

ERL_NIF_TERM
erl_nmatrix_from_list(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (2!=argc) {
        return enif_make_badarg(env);
    }

    const char * error = 0;
    create_options_t options;
    if (!parse_create_options(env, argv[1], &options, 0, &error)) {
        return make_error(env, error);
    }

    if (!enif_is_list(env, argv[0])) {
        return make_error(env, ERR_ARG_BAD_LIST);
    }

    ERL_NIF_TERM row_hd = 0, row_tl = argv[0];
    ERL_NIF_TERM col_hd = 0, col_tl = 0;
    unsigned rows_count = 0;
    enif_get_list_length(env, argv[0], &rows_count);
    unsigned max_cols_count = 0;
    unsigned cols_count = 0;

    // Calculate required space
    while (enif_get_list_cell(env, row_tl, &row_hd, &row_tl)) {
        if (!enif_is_list(env, row_hd)) {
            return make_error(env, ERR_ARG_BAD_2D_LIST);
        }
        enif_get_list_length(env, row_hd, &cols_count);
        max_cols_count = cols_count > max_cols_count ? cols_count : max_cols_count;
    }

    const size_t elements_count = rows_count * max_cols_count;
    scalar_element_t * elements = enif_alloc(elements_count*sizeof(*elements));
    memset(elements, 0, elements_count*sizeof(*elements));
    scalar_element_t * current = 0;
    unsigned element_index = 0;
    unsigned row_index = 0;
    unsigned col_index = 0;
    _Bool error_occured = 0;
    row_tl = argv[0];
    while (enif_get_list_cell(env, row_tl, &row_hd, &row_tl)) {
        col_tl = row_hd;
        col_index = 0;
        col_hd = 0;
        while (enif_get_list_cell(env, col_tl, &col_hd, &col_tl)) {
            element_index = row_index * max_cols_count + col_index;
            current = elements + element_index*sizeof(*current);
            if (!parse_erl_scalar_term(env, col_hd, current, &options, &error)) {
                error_occured = 1;
                break;
            }
            col_index++;
        }
        while (col_index < max_cols_count) {
            element_index = row_index * max_cols_count + col_index;
            current = elements + element_index*sizeof(*current);
            current->dtype = DTSingle;
            current->value.s = 0.0;
            col_index++;
        }
        if (error_occured)
            break;
        row_index++;
    }
    if (error_occured) {
        if (elements)
            enif_free(elements);
        return make_error(env, error);
    }

    nmatrix_t result;
    memset(&result, 0, sizeof(result));
    result.view0.offset = result.view1.offset = 0;
    result.view0.increment = (int32_t) max_cols_count;
    result.view1.increment = 1;
    result.view0.size = result.array.size0 = rows_count;
    result.view1.size = result.array.size1 = max_cols_count;
    result.array.dtype = common_dtype_for_elements(elements_count, elements);

    const size_t item_size = scalar_size(result.array.dtype);
    const size_t bin_size = item_size * result.array.size0 * result.array.size1;
    if (!enif_alloc_binary(bin_size, &result.array.bin)) {
        return make_error(env, ERR_GEN_MEMORY_ALLOCATION);
    }
    void * p = result.array.bin.data;
    for (size_t i=0; i<elements_count; ++i) {
        current = elements + i*sizeof(*current);
        scalar_value_t fill_value = convert_type(current->value, current->dtype, result.array.dtype);
        memcpy(p, &fill_value, item_size);
        p += item_size;
    }
    enif_free(elements);
    return nmatrix_to_erl_record(env, &result, 0);
}

ERL_NIF_TERM
erl_nmatrix_row(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (2!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nmatrix_t mat;
    ERL_NIF_TERM bin_ref = 0;
    if (!parse_erl_matrix(env, argv[0], &mat, &error, &bin_ref)) {
        return make_error(env, error);
    }
    unsigned row_no = 0;
    if (!enif_get_uint(env, argv[1], &row_no)) {
        return make_error(env, ERR_ARG_BAD_INDEX);
    }
    if (row_no >= mat.view0.size) {
        return make_error(env, ERR_ARG_INDEX_OUT_OF_RANGE);
    }
    nvector_t result;
    result.array.size0 = mat.array.size0;
    result.array.size1 = mat.array.size1;
    memset(&result.array.bin, 0, sizeof(result.array.bin));
    result.array.dtype = mat.array.dtype;
    result.view.offset =
            mat.view0.offset +
            mat.view0.increment*row_no;
    result.view.size = mat.view1.size;
    result.view.increment = mat.view1.increment;

    return nvector_to_erl_record(env,
                                 &result,
                                 enif_make_copy(env, bin_ref)
                                 );
}

ERL_NIF_TERM
erl_nmatrix_col(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (2!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nmatrix_t mat;
    ERL_NIF_TERM bin_ref = 0;
    if (!parse_erl_matrix(env, argv[0], &mat, &error, &bin_ref)) {
        return make_error(env, error);
    }
    unsigned col_no = 0;
    if (!enif_get_uint(env, argv[1], &col_no)) {
        return make_error(env, ERR_ARG_BAD_INDEX);
    }
    if (col_no >= mat.view1.size) {
        return make_error(env, ERR_ARG_INDEX_OUT_OF_RANGE);
    }
    nvector_t result;
    result.array.size0 = mat.array.size0;
    result.array.size1 = mat.array.size1;
    memset(&result.array.bin, 0, sizeof(result.array.bin));
    result.array.dtype = mat.array.dtype;
    result.view.offset =
            mat.view1.offset +
            mat.view1.increment*col_no;
    result.view.size = mat.view0.size;
    result.view.increment = mat.view0.increment;

    return nvector_to_erl_record(env,
                                 &result,
                                 enif_make_copy(env, bin_ref)
                                 );
}

ERL_NIF_TERM
erl_nmatrix_diag(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (2!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nmatrix_t mat;
    ERL_NIF_TERM bin_ref = 0;
    if (!parse_erl_matrix(env, argv[0], &mat, &error, &bin_ref)) {
        return make_error(env, error);
    }
    int diag_no = 0;
    if (!enif_get_int(env, argv[1], &diag_no)) {
        return make_error(env, ERR_ARG_BAD_DIAGONAL);
    }
    nvector_t result;
    result.array.size0 = mat.array.size0;
    result.array.size1 = mat.array.size1;
    memset(&result.array.bin, 0, sizeof(result.array.bin));
    result.array.dtype = mat.array.dtype;
    unsigned udiag_no = 0;
    if (diag_no >= 0) {
        udiag_no = diag_no;
        result.view.size = udiag_no < mat.view1.size
                ? mat.view1.size - udiag_no : 0;
        result.view.size = min_size(result.view.size, mat.view0.size);
        result.view.offset = mat.view1.increment*udiag_no;
    }
    else {
        udiag_no = -diag_no;
        result.view.size = udiag_no < mat.view0.size
                ? mat.view0.size - udiag_no : 0;
        result.view.size = min_size(result.view.size, mat.view1.size);
        result.view.offset = mat.view0.increment*udiag_no;
    }
    result.view.offset += mat.view0.offset + mat.view1.offset;
    result.view.increment = mat.view0.increment + mat.view1.increment;

    return nvector_to_erl_record(env,
                                 &result,
                                 enif_make_copy(env, bin_ref)
                                 );
}

ERL_NIF_TERM
erl_nmatrix_get(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (3!=argc) {
        return enif_make_badarg(env);
    }
    const char* error = 0;
    nmatrix_t mat;
    ERL_NIF_TERM bin_ref = 0;
    if (!parse_erl_matrix(env, argv[0], &mat, &error, &bin_ref)) {
        return make_error(env, error);
    }

    size_t indeces[2];
    if (!parse_erl_shape(env, argv[1], 2, indeces, &error)) {
        return make_error(env, ERR_ARG_BAD_INDEX);
    }
    if (indeces[0] >= mat.view0.size || indeces[1] >= mat.view1.size) {
        return make_error(env, ERR_ARG_INDEX_OUT_OF_RANGE);
    }
    convert_option_t convert;
    if (!parse_convert_option(env, argv[2], &convert, &error)) {
        return make_error(env, error);
    }

    scalar_element_t element;
    element.dtype = mat.array.dtype;
    element.value = nmatrix_get(&mat, indeces[0], indeces[1]);
    ERL_NIF_TERM result = make_scalar_value(env, element, convert);

    return result;
}

scalar_value_t
nmatrix_get(const nmatrix_t *mat, size_t row_index, size_t col_index)
{
    scalar_value_t element;
    const size_t item_size = scalar_size(mat->array.dtype);
    const size_t row_offset =
            mat->view0.offset + row_index * mat->view0.increment;
    const size_t col_offset =
            mat->view1.offset + col_index * mat->view1.increment;
    const size_t real_index = row_offset + col_offset;
    void * src = mat->array.bin.data + item_size * real_index;
    memcpy(&element, src, item_size);
    return element;
}
