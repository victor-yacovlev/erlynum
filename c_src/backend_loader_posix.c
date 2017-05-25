#include "backend_loader.h"
#include "error_atoms.h"
#ifdef _POSIX_SOURCE

#include <errno.h>
#include <dlfcn.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

static void * Lib_Ptr = 0;

static _Bool
take_next_path(const char * s, char *head, size_t head_size, const char* *tail)
{
    memset(head, 0, head_size);
    if (!s) {
        return false;
    }
    const char * colon = 0;
    if ((colon=strstr(s, ":"))) {
        if (colon!=s) {
            strncpy(head, s, colon-s);
        }
        *tail = colon + 1;
        return true;
    }
    else {
        if (strlen(s) > 0) {
            strncpy(head, s, strlen(s));
            *tail = 0;
            return true;
        }
    }
    return false;
}

static _Bool
find_library_in_path(const char* names[], const char *canonical_name, char* result)
{
    char complete_path_env[4096];
    memset(complete_path_env, 0, sizeof(complete_path_env));
    const char * path_env = getenv("LD_LIBRARY_PATH");
    size_t add_pos = 0;
    if (path_env) {
        strncat(complete_path_env, path_env, sizeof(complete_path_env)-add_pos);
        add_pos = strlen(complete_path_env);
    }
    if (8==sizeof(void*)) {
        strncat(complete_path_env, ":/lib64", sizeof(complete_path_env)-add_pos);
        add_pos = strlen(complete_path_env);
        strncat(complete_path_env, ":/usr/lib64", sizeof(complete_path_env)-add_pos);
        add_pos = strlen(complete_path_env);
        strncat(complete_path_env, ":/usr/local/lib64", sizeof(complete_path_env)-add_pos);
        add_pos = strlen(complete_path_env);
    }
    else {
        strncat(complete_path_env, ":/lib", sizeof(complete_path_env)-add_pos);
        add_pos = strlen(complete_path_env);
        strncat(complete_path_env, ":/usr/lib", sizeof(complete_path_env)-add_pos);
        add_pos = strlen(complete_path_env);
        strncat(complete_path_env, ":/usr/local/lib", sizeof(complete_path_env)-add_pos);
        add_pos = strlen(complete_path_env);
    }
    const char * tail = complete_path_env;
    char prefix[PATH_MAX];
    char buffer[PATH_MAX];
    struct stat st;
    size_t pref_len = 0;
    while (take_next_path(tail, prefix, sizeof(prefix), &tail)) {
        pref_len = strlen(prefix);
        if (pref_len > 0) {
            for (const char* *it = names; *it!=NULL; ++it) {
                const char * name = *it;
                memset(buffer, 0, sizeof(buffer));
                strncpy(buffer, prefix, sizeof(buffer));
                buffer[pref_len] = '/';
                strncpy(buffer + pref_len + 1, name, strlen(name));
                memset(&st, 0, sizeof(st));
                if (0==stat(buffer, &st)) {
                    if (S_ISREG(st.st_mode) || S_ISLNK(st.st_mode)) {
                        strncpy(result, buffer, PATH_MAX);
                        return true;
                    }
                }

                memset(buffer, 0, sizeof(buffer));
                strncpy(buffer, prefix, sizeof(buffer));
                const size_t can_size = strlen(canonical_name);
                buffer[pref_len] = '/';
                strncpy(buffer + pref_len + 1, canonical_name, can_size);
                buffer[pref_len+can_size+1] = '/';
                strncpy(buffer+pref_len+can_size+2, name, strlen(name));
                memset(&st, 0, sizeof(st));
                if (0==stat(buffer, &st)) {
                    if (S_ISREG(st.st_mode) || S_ISLNK(st.st_mode)) {
                        strncpy(result, buffer, PATH_MAX);
                        return true;
                    }
                }
            }
        }
    }
    return false;
}


static _Bool
load_library(const char * file_name, _Bool global_symbols,
             void* *out,
             const char* *error_reason)
{
    const int mode = global_symbols
            ? RTLD_GLOBAL | RTLD_NOW
            : RTLD_NOW;
    *out = dlopen(file_name, mode);
    if (*out) {
        return true;
    }
    else {
        *error_reason = dlerror();
        return false;
    }
}

static void
unload_library(void * lib)
{
    if (lib) {
        dlclose(lib);
    }
}

static _Bool
load_atlas(const char* *error, const char* *error_reason)
{
    static const char * Names[] = {
        "libsatlas.so.3", "libsatlas.so", "libatlas.so.3", "libatlas.so", NULL
    };
    char atlas_path[PATH_MAX];
    memset(atlas_path, 0, sizeof(atlas_path));
    void * atlas_lib = 0;

    if (!find_library_in_path(Names, "atlas", atlas_path)) {
        *error = ERR_GEN_NOTFOUND_ATLAS;
        return false;
    }

    if (!load_library(atlas_path, false, &atlas_lib, error_reason)) {
        *error = ERR_GEN_DLOPEN;
        return false;
    }

    Lib_Ptr = atlas_lib;
    return true;
}

static _Bool
load_blas(const char* *error, const char* *error_reason)
{
    static const char * Blas_Names[] = {
        "libblas.so.3", "libblas.so", NULL
    };

    static const char * CBlas_Names[] = {
        "libcblas.so.3", "libcblas.so", NULL
    };
    *error = 0;
    *error_reason = 0;
    char blas_path[PATH_MAX];
    memset(blas_path, 0, sizeof(blas_path));
    char cblas_path[PATH_MAX];
    memset(cblas_path, 0, sizeof(blas_path));
    void * blas_lib = 0;
    void * cblas_lib = 0;

    if (!find_library_in_path(Blas_Names, "blas", blas_path)) {
        *error = ERR_GEN_NOTFOUND_BLAS;
        return false;
    }

    if (!find_library_in_path(CBlas_Names, "blas", cblas_path)) {
        *error = ERR_GEN_NOTFOUND_CBLAS;
        return false;
    }

    if (!load_library(blas_path, true, &blas_lib, error_reason)) {
        *error = ERR_GEN_DLOPEN;
        return false;
    }

    if (!load_library(cblas_path, false, &cblas_lib, error_reason)) {
        unload_library(blas_lib);
        *error = ERR_GEN_DLOPEN;
        return false;
    }

    Lib_Ptr = cblas_lib;
    return true;
}

static _Bool
load_intel_mkl(const char* *error, const char* *error_reason)
{
    static const char * MklRT_Names[] = { "libmkl_rt.so", NULL };
    char mkl_rt_path[PATH_MAX];
    void * mkl_rt_lib = 0;
    if (!find_library_in_path(MklRT_Names, "mkl", mkl_rt_path)) {
        *error = ERR_GEN_NOTFOUND_BLAS;
    }
    if (!load_library(mkl_rt_path, false, &mkl_rt_lib, error_reason)) {
        *error = ERR_GEN_DLOPEN;
        return false;
    }
    Lib_Ptr = mkl_rt_lib;
    return true;
}

static void
prepend_to_ld_library_path(const char *path)
{
    char * old_path = getenv("LD_LIBRARY_PATH");
    static char new_path[512];
    memset(new_path, 0, sizeof(new_path));
    const size_t priv_len = strlen(path);
    memcpy(new_path, path, priv_len);
    if (old_path) {
        new_path[priv_len] = ':';
        strncpy(new_path+priv_len+1, old_path, sizeof(new_path)-1-priv_len);
    }
    setenv("LD_LIBRARY_PATH", new_path, 1);
}

_Bool
init_backend(const char *priv_dir, const char *backend_name,
             const char* *error, const char* *error_reason)
{
    if (8==sizeof(void*)) {
        prepend_to_ld_library_path("/opt/intel/mkl/lib/intel64");
    }
    else {
        prepend_to_ld_library_path("/opt/intel/mkl/lib/ia32");
    }
    prepend_to_ld_library_path(priv_dir);
    static char local_error[255];
    memset(local_error, 0, sizeof(local_error));
    strncpy(local_error, backend_name, sizeof(local_error));
    if (0==strncmp("atlas", backend_name, 5)) {
        return load_atlas(error, error_reason);
    }
    else if (0==strncmp("blas", backend_name, 4)) {
        return load_blas(error, error_reason);
    }
    else if (0==strncmp("mkl", backend_name, 4)) {
        return load_intel_mkl(error, error_reason);
    }
    else if (0==strncmp("auto", backend_name, 4)) {
        if (load_intel_mkl(error, error_reason))    return true;
        if (load_atlas(error, error_reason))        return true;
        if (load_blas(error, error_reason))         return true;
    }
    *error_reason = local_error;
    *error = ERR_GEN_BACKEND_LOAD;
    return false;
}

void*
resolve_blas_function(const char* func_name)
{
    if (!Lib_Ptr)
        return 0;    
    return dlsym(Lib_Ptr, func_name);
}

#endif // _POSIX_SOURCE
