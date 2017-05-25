%% @author          Victor Yacovlev <v.yacovlev@gmail.com>
%% @copyright       2017 Victor Yacovlev

%% @doc Matrix operation functions.
%%
%% Matrix is a 2-dimensional array, which elements stored sequentially into Erlang binary object of
%% some IEEE758 floating point data type: single float, double float, single complex or double complex.
%%
%% Note that matrix indexes starts from 0, but not from 1 as for lists.

-module(nmatrix).
-include("erlynum.hrl").
-define(WE, erlynum_p:wrap_error).
-export([zeros/1, zeros/2, ones/1, ones/2, full/2, full/3]).
-export([identity/1, identity/2, eye/1, eye/2, eye/3, eye/4]).
-export([to_list/1, to_list/2, from_list/1, from_list/2]).
-export([get/2, get/3, row/2, col/2, diag/1, diag/2]).
-export([transpose/1]).

-spec zeros(Shape :: erlynum:shape()) -> erlynum:nmatrix().
%% @equiv zeros(Shape, [])
zeros(Shape) -> zeros(Shape, []).

-spec zeros(
    Shape           ::  erlynum:shape(),
    Options         ::  [ erlynum:create_option() ]
) -> erlynum:nmatrix().
%% @doc Returns a matrix of the given `Shape' filled by zeroes.
zeros(Shape, Options) -> ?WE(erlynum_nif:nmatrix_full(Shape, 0, Options)).

-spec ones(Shape :: erlynum:shape()) -> erlynum:nmatrix().
%% @equiv ones(Shape, [])
ones(Shape) -> ones(Shape, []).

-spec ones(
    Shape           ::  erlynum:shape(),
    Options         ::  [ erlynum:create_option() ]
) -> erlynum:nmatrix().
%% @doc Returns a matrix of the given `Shape' filled by ones.
ones(Shape, Options) -> ?WE(erlynum_nif:nmatrix_full(Shape, 1, Options)).

-spec full(
    Shape           ::  erlynum:shape(),
    FillValue       ::  erlynum:nscalar()
) -> erlynum:nmatrix().
%% @equiv full(Shape, FillValue, [])
full(Shape, FillValue) -> ?WE(erlynum_nif:nmatrix_full(Shape, FillValue, [])).

-spec full(
    Shape           ::  erlynum:shape(),
    FillValue       ::  erlynum:nscalar(),
    Options         ::  [ erlynum:create_option() ]
) -> erlynum:nmatrix().
%% @doc Returns a matrix of the given `Shape' filled by the specified scalar value `FillValue'.
full(Shape, FillValue, Options) -> ?WE(erlynum_nif:nmatrix_full(Shape, FillValue, Options)).

-spec identity(N :: non_neg_integer()) -> erlynum:nmatrix().
%% @equiv identity(N, [])
identity(N) -> identity(N, []).

-spec identity(
    N               ::  non_neg_integer(),
    Options         ::  [ erlynum:create_option() ]
) -> erlynum:nmatrix().
%% @doc Returns a square matrix of size `N' with ones on the main diagonal.
identity(N, Options) -> eye({N, N}, 1.0, 0, Options).

-spec eye(Shape :: erlynum:shape()) -> erlynum:nmatrix().
%% @doc Returns a matrix of the given `Shape' filled by ones on
%% the main diagonal and zeroes in all other places.
eye(Shape) -> eye(Shape, 1.0).

-spec eye(
    Shape           ::  erlynum:shape(),
    DiagValue_Opts  ::  erlynum:nscalar()
                    |   [ erlynum:create_option() ]
) -> erlynum:nmatrix().
%% @doc Returns a matrix of the given `Shape' filled by ones of specified value on
%% the main diagonal and zeroes in all other places.
%%
%% If last argument is a scalar value, it specifies a value to be filled.
%% If last argument is a list, it specifies create options.
eye(Shape, DiagValue_Opts) when is_list(DiagValue_Opts) ->
    eye(Shape, 1.0, 0, DiagValue_Opts);
eye(Shape, DiagValue_Opts) ->
    eye(Shape, DiagValue_Opts, 0).

-spec eye(
    Shape           ::  erlynum:shape(),
    DiagonalValue   ::  erlynum:nscalar(),
    DiagNum_Options ::  integer()
    |   [ erlynum:create_option() ]
) -> erlynum:nmatrix().
%% @doc Returns a matrix of the given `Shape' filled by specified scalar `DiagonalValue' on
%% the specified of main diagonal and zeroes in all other places.
%%
%% If last argument is an integer value, it specifies diagonal number.
%% If last argument is a list, it specifies create options.
eye(Shape, DiagonalValue, DiagNum_Options) when is_integer(DiagNum_Options) ->
    eye(Shape, DiagonalValue, DiagNum_Options, []);
eye(Shape, DiagonalValue, DiagNum_Options) when is_list(DiagNum_Options) ->
    eye(Shape, DiagonalValue, 0, DiagNum_Options).

-spec eye(
    Shape           ::  erlynum:shape(),
    DiagonalValue   ::  erlynum:nscalar(),
    DiagonalNumber  ::  integer(),
    Options         ::  [ erlynum:create_option() ]
) -> erlynum:nmatrix().
%% @doc Returns a matrix of the given `Shape' filled by specified scalar `DiagonalValue' on
%% the specified diagonal `DiagonalNumber' and zeroes in all other places.
eye(Shape, DiagonalValue, DiagonalNumber, Options) ->
    ?WE(erlynum_nif:nmatrix_eye(Shape, DiagonalValue, DiagonalNumber, Options)).

-spec from_list([[erlynum:nscalar()]]) -> erlynum:nmatrix().
%% @equiv from_list(List, [])
from_list(List) -> from_list(List, []).

-spec from_list([[erlynum:nscalar()]], [erlynum:create_option()]) -> erlynum:nmatrix().
%% @doc Returns a matrix from the given two-dimensional `List' items.
from_list(List, Options) -> ?WE(erlynum_nif:nmatrix_from_list(List, Options)).

-spec to_list(NMatrix :: erlynum:nmatrix()) ->  [ [ erlynum:nscalar() ] ].
%% @equiv to_list(NMatrix, [])
to_list(NMatrix) -> to_list(NMatrix, noconvert).

-spec to_list(
    NMatrix         ::  erlynum:nmatrix(),
    Convert         ::  erlynum:convert_option()
) ->  [ [ erlynum:nscalar() ] ].
%% @doc Returns an Erlang list from the given matrix.
to_list(NMatrix, Convert) -> ?WE(erlynum_nif:nmatrix_to_list(NMatrix, Convert)).

-spec get(
    NMatrix         ::  erlynum:nmatrix(),
    Index           ::  { non_neg_integer(), non_neg_integer() }
) -> erlynum:nscalar().
%% @equiv get(NMatrix, Index, noconvert)
get(NMatrix, Index) -> get(NMatrix, Index, noconvert).

-spec get(
    NMatrix         ::  erlynum:nmatrix(),
    Index           ::  { non_neg_integer(), non_neg_integer() },
    Convert         ::  erlynum:convert_option()
) -> erlynum:nscalar().
%% @doc Returns a scalar element value from the matrix.
%% Note that `Index' starts from 0, but not 1.
get(NMatrix, Index, Convert) -> ?WE(erlynum_nif:nmatrix_get(NMatrix, Index, Convert)).

-spec row(
    NMatrix         ::  erlynum:nmatrix(),
    RowIndex        ::  non_neg_integer()
) -> erlynum:nvector().
%% @doc Returns a row from the matrix as vector view.
%% Note that `RowIndex' starts from 0, but not 1.
row(NMatrix, RowIndex) -> ?WE(erlynum_nif:nmatrix_row(NMatrix, RowIndex)).

-spec col(
    NMatrix         ::  erlynum:nmatrix(),
    ColIndex        ::  non_neg_integer()
) -> erlynum:nvector().
%% @doc Returns a column from the matrix as vector view.
%% Note that `ColIndex' starts from 0, but not 1.
col(NMatrix, ColIndex) -> ?WE(erlynum_nif:nmatrix_col(NMatrix, ColIndex)).

-spec diag(NMatrix :: erlynum:nmatrix()) -> erlynum:nvector().
%% @doc Returns the main diagonal from the matrix as vector view.
diag(NMatrix) -> diag(NMatrix, 0).

-spec diag(
    NMatrix         ::  erlynum:nmatrix(),
    DiagNumber      ::  integer()
) -> erlynum:nvector().
%% @doc Returns a specified diagonal from the matrix as vector view.
%% Main diagonal number is 0. `DiagNumber' may be both positive and negative.
diag(NMatrix, DiagNumber) -> ?WE(erlynum_nif:nmatrix_diag(NMatrix, DiagNumber)).

-spec transpose(NMatrix :: erlynum:nmatrix()) -> erlynum:nmatrix().
%% @doc Returns a transposed view of matrix.
transpose(NMatrix) ->
    {RowView, ColView} = NMatrix#nmatrix.view,
    NMatrix#nmatrix{view = {ColView, RowView}}.

