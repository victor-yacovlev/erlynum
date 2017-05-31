#include "backend_loader.h"
#include "nvector.h"
#include "nmatrix.h"
#include "nblas.h"


#include <erl_nif.h>

static ErlNifFunc
nif_funcs[] = {
    {"init_backend", 2, erl_init_backend},
    {"backend_in_use", 0, erl_backend_in_use},

    {"nvector_full", 3, erl_nvector_full},
    {"nvector_from_list", 2, erl_nvector_from_list},
    {"nvector_to_list", 2, erl_nvector_to_list},
    {"nvector_get", 3, erl_nvector_get},
    {"nvector_range", 4, erl_nvector_range},
    {"nvector_linspace", 4, erl_nvector_linspace},
    {"nvector_logspace", 5, erl_nvector_logspace},
    {"nvector_geomspace", 4, erl_nvector_geomspace},
    {"nvector_copy", 2, erl_nvector_copy},
    {"nvector_scale", 3, erl_nvector_scale},
    {"nvector_axpy", 4, erl_nvector_axpy},
    {"nvector_dot", 3, erl_nvector_dot},
    {"nvector_asum", 2, erl_nvector_asum},
    {"nvector_iamax_iamin", 2, erl_nvector_iamax_iamin},
    {"nvector_nrm2", 2, erl_nvector_nrm2},

    {"nmatrix_full", 3, erl_nmatrix_full},
    {"nmatrix_to_list", 2, erl_nmatrix_to_list},
    {"nmatrix_from_list", 2, erl_nmatrix_from_list},
    {"nmatrix_eye", 4, erl_nmatrix_eye},

    {"nmatrix_row", 2, erl_nmatrix_row},
    {"nmatrix_col", 2, erl_nmatrix_col},
    {"nmatrix_diag", 2, erl_nmatrix_diag},
    {"nmatrix_get", 3, erl_nmatrix_get}

};

ERL_NIF_INIT(erlynum_nif, nif_funcs, NULL, NULL, NULL, NULL)
