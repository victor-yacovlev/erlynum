%% @author          Victor Yacovlev <v.yacovlev@gmail.com>
%% @copyright       2017 Victor Yacovlev
%% @private
-module(erlynum_nif).
-define(APPNAME, erlynum).

-on_load(init/0).

-export([
    nvector_full/3, nvector_from_list/2, nvector_to_list/2,
    nvector_get/3, nvector_range/4, nvector_linspace/4,
    nvector_logspace/5, nvector_geomspace/4,
    nmatrix_full/3, nmatrix_to_list/2, nmatrix_from_list/2,
    nmatrix_get/3, nmatrix_eye/4, nmatrix_row/2, nmatrix_col/2,
    nmatrix_diag/2, nvector_copy/2, nvector_scale/3, nvector_axpy/4, nvector_dot/3,
    nvector_asum/2, nvector_iamax_iamin/2]).


init() ->
    Backend = application:get_env(?APPNAME, backend, auto),
    case load_driver() of
        ok                      ->
            case init_backend(code:priv_dir(?APPNAME), Backend) of
                ok              -> ok;
                Other           -> Other
            end;
        {error, {upgrade, _}}   -> ok;
        {error, {Reason, Text}} ->
            io:format("Error loading erlynum driver: ~p, ~s~n", [Reason, Text]),
            {error, Reason}
    end.


load_driver() ->
    PrivDir = code:priv_dir(?APPNAME),
    DriverName = filename:join(PrivDir, atom_to_list(?APPNAME) ++ "_drv"),
    erlang:load_nif(DriverName, 0).

-define(NL, {error, erlynum_drv_not_loaded}).

init_backend(_PrivDir, _Backend) -> ?NL.

nvector_full(_Shape, _FillValue, _Options) -> ?NL.
nvector_from_list(_Contents, _Options) -> ?NL.
nvector_to_list(_Vector, _Convert) -> ?NL.
nvector_get(_Vector, _Index, _Convert) -> ?NL.
nvector_range(_Start, _Stop, _Step, _Options) -> ?NL.
nvector_linspace(_Start, _Stop, _Count, _Options) -> ?NL.
nvector_logspace(_Start, _Stop, _Count, _Base, _Options) -> ?NL.
nvector_geomspace(_Start, _Stop, _Count, _Options) -> ?NL.
nvector_copy(_Vector, _Options) -> ?NL.
nvector_scale(_Vector, _Alpha, _Options) -> ?NL.
nvector_axpy(_Y, _X, _Alpha, _Options) -> ?NL.
nvector_dot(_Y, _X, _Options) -> ?NL.
nvector_asum(_Vector, _Options) -> ?NL.
nvector_iamax_iamin(_Vector, _Mode) -> ?NL.

nmatrix_full(_Shape, _FillValue, _Options) -> ?NL.
nmatrix_to_list(_Matrix, _Convert) -> ?NL.
nmatrix_from_list(_Contents, _Options) -> ?NL.
nmatrix_eye(_Shape, _DiagValue, _DiagNumber, _Options) -> ?NL.
nmatrix_get(_Matrix, _Index, _Convert) -> ?NL.
nmatrix_row(_Matrix, _RowNumber) -> ?NL.
nmatrix_col(_Matrix, _ColNumber) -> ?NL.
nmatrix_diag(_NMatrix, _DiagNumber) -> ?NL.