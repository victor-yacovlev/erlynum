#ifndef ERROR_ATOMS_H
#define ERROR_ATOMS_H

// Generic errors ERR_GEN_...
#define ERR_GEN_MEMORY_ALLOCATION       "memory_allocation_error"
#define ERR_GEN_NOTFOUND_ATLAS          "libatlas_not_found"
#define ERR_GEN_NOTFOUND_BLAS           "libblas_not_found"
#define ERR_GEN_NOTFOUND_CBLAS          "libcblas_not_found"
#define ERR_GEN_NOTFOUND_MKL_RT         "libmkl_rt_not_found"
#define ERR_GEN_DLOPEN                  "dlopen_error"
#define ERR_GEN_BACKEND_LOAD            "backend_load_error"

// BLAS backend errors ERR_BLAS_...
#define ERR_BLAS_NOTFOUND_COPY          "cblas_copy_not_resolved_in_backend"
#define ERR_BLAS_NOTFOUND_SCAL          "cblas_scal_not_resolved_in_backend"
#define ERR_BLAS_NOTFOUND_AXPY          "cblas_axpy_not_resolved_in_backend"
#define ERR_BLAS_NOTFOUND_DOT           "cblas_dot_not_resolved_in_backend"
#define ERR_BLAS_NOTFOUND_SDOT          "cblas_sdot_not_resolved_in_backend"
#define ERR_BLAS_NOTFOUND_ASUM          "cblas_asum_not_resolved_in_backend"
#define ERR_BLAS_NOTFOUND_IAMIN         "cblas_iamin_not_resolved_in_backend"
#define ERR_BLAS_NOTFOUND_IAMAX         "cblas_iamax_not_resolved_in_backend"
#define ERR_BLAS_NOTSUPPORT_DTYPE       "cblas_unsupported_dtype"

// Errors while parsing broken records or erlynum structures ERR_REC_...
#define ERR_REC_VIEW                    "bad_view_parameters"
#define ERR_REC_NVECTOR                 "bad_nvector"
#define ERR_REC_NMATRIX                 "bad_nmatrix"
#define ERR_REC_DTYPE                   "bad_dtype_in_record"
#define ERR_REC_BIN                     "bad_binary_in_record"

// Wrong argument value errors ERR_ARG_...
#define ERR_ARG_BAD_INDEX               "bad_index_value"
#define ERR_ARG_BAD_DIAGONAL            "bad_diagonal_value"
#define ERR_ARG_BAD_LIST                "bad_list_value"
#define ERR_ARG_BAD_2D_LIST             "bad_2dimension_list_value"
#define ERR_ARG_INDEX_OUT_OF_RANGE      "index_out_of_range"
#define ERR_ARG_BAD_COUNT               "bad_count_value"
#define ERR_ARG_BAD_START               "bad_start_value"
#define ERR_ARG_BAD_STOP                "bad_stop_value"
#define ERR_ARG_BAD_STEP                "bad_step_value"
#define ERR_ARG_BAD_RANGE               "bad_range_value"
#define ERR_ARG_BAD_NUMBER              "bad_number_value"
#define ERR_ARG_BAD_BOOL                "bad_boolean_value"
#define ERR_ARG_BAD_NUMBER_OR_COMPLEX   "bad_number_or_complex_value"
#define ERR_ARG_BAD_PROPLIST            "bad_property_list"
#define ERR_ARG_BAD_CREATE_OPTION       "bad_create_option"
#define ERR_ARG_BAD_DTYPE               "bad_dtype_value"
#define ERR_ARG_BAD_PRECISION           "bad_precision_value"
#define ERR_ARG_BAD_SHAPE               "bad_shape_value"
#define ERR_ARG_BAD_CONVERT             "bad_convert_value"
#define ERR_ARG_NOT_REAL_ARRAY          "not_real_value_array"
#define ERR_ARG_NOT_REAL_SCALAR         "not_real_value_scalar"
#define ERR_ARG_BAD_MIN_MAX_ATOM        "bad_min_or_max_mode_atom"

#endif // ERROR_ATOMS_H

