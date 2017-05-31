#include "backend_loader.h"
#include "erl_util.h"

#include <string.h>

ERL_NIF_TERM
erl_init_backend(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    static const char * DefaultBackendName = "cblas";
    const char * error = 0;
    const char * error_reason = 0;
    char backend_name[16];      memset(backend_name, 0, sizeof(backend_name));
    char priv_dir[256];         memset(priv_dir, 0, sizeof(priv_dir));
    if (1==argc) {
        strncpy(backend_name, DefaultBackendName, sizeof(backend_name));
    }
    else if (2==argc && enif_is_atom(env, argv[1])) {
        enif_get_atom(env, argv[1], backend_name, sizeof(backend_name), ERL_NIF_LATIN1);
    }
    else {
        return enif_make_badarg(env);
    }
    if (!enif_get_string(env, argv[0], priv_dir, sizeof(priv_dir), ERL_NIF_LATIN1)) {
        return enif_make_badarg(env);
    }
    if (init_backend(priv_dir, backend_name, &error, &error_reason)) {
        return make_atom(env, "ok");
    }
    else {
        return error_reason
                ? make_error_with_reason(env, error, error_reason)
                : make_error(env, error);
    }
}

ERL_NIF_TERM
erl_backend_in_use(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[])
{
    if (0!=argc) {
        return enif_make_badarg(env);
    }
    const blas_backend_driver_t driver = backend_in_use();
    const char * name = 0;
    switch (driver) {
    case DrvNone:       name = "no_backend";    break;
    case DrvIntelMKL:   name = "mkl";           break;
    case DrvATLAS:      name = "atlas";         break;
    case DrvNetlibBLAS: name = "blas";          break;
    default:            name = "other";         break;
    }
    return make_atom(env, name);
}
