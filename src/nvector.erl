%% @author          Victor Yacovlev <v.yacovlev@gmail.com>
%% @copyright       2017 Victor Yacovlev

%% @doc Vector operation functions.
%%
%% Vector is a 1-dimensional array, which elements stored sequentially into Erlang binary object of
%% some IEEE758 floating point data type: single float, double float, single complex or double complex.
%%
%% Vector contents uses regular Erlang binary storage and relies on it's garbage collection. In some
%% cases it is possible to refer the same binary contents.
%%
%% Note that vector indexes starts from 0, but not from 1 as for lists.
%% @end

-module(nvector).
-include_lib("eunit/include/eunit.hrl").
-include("erlynum.hrl").

-export([full/2, full/3, zeros/1, zeros/2, ones/1, ones/2]).
-export([range/1, range/2, range/3, range/4]).
-export([linspace/3, linspace/4, logspace/3, logspace/4, logspace/5, geomspace/3, geomspace/4]).
-export([from_list/1, from_list/2, to_list/1, to_list/2]).
-export([get/2, get/3]).
-export([
    copy/1, copy/2,
    sum/2, sum/3,
    scale/2, scale/3,
    axpy/3, axpy/4,
    dot/3, dot/2,
    asum/2, asum/1,
    iamax/1, iamin/1
    , nrm2/2, nrm2/1, euclidean_norm/1, euclidean_norm/2]).

-define(WE, erlynum_p:wrap_error).

-spec zeros(non_neg_integer()) -> erlynum:nvector().
%% @equiv zeros(Size, [{dtype, d}])
zeros(Size) -> zeros(Size, [{dtype, d}]).

-spec zeros(non_neg_integer(), [erlynum:create_option()]) -> erlynum:nvector().
%% @doc Returns a vector of the given size filled by zeroes.
zeros(Size, Options) -> full(Size, 0, Options).

-spec ones(non_neg_integer()) -> erlynum:nvector().
%% @equiv ones(Size, [{dtype, d}])
ones(Size) -> ones(Size, [{dtype, d}]).

-spec ones(non_neg_integer(), [erlynum:create_option()]) -> erlynum:nvector().
%% @doc Returns a vector of the given size filled by ones.
ones(Size, Options) -> ?WE(erlynum_nif:nvector_full(Size, 1, Options)).

-spec full(non_neg_integer(), erlynum:nscalar()) -> erlynum:nvector().
%% @equiv full(Size, InitialValue, [{dtype, auto}])
full(Size, InitialValue) -> ?WE(erlynum_nif:nvector_full(Size, InitialValue, [{dtype, auto}])).

-spec full(non_neg_integer(), erlynum:nscalar(), [erlynum:create_option()]) -> erlynum:nvector().
%% @doc Returns a vector of the given size filled by the specified scalar value `InitialValue'.
full(Size, InitialValue, Options) -> ?WE(erlynum_nif:nvector_full(Size, InitialValue, Options)).

-spec from_list([erlynum:nscalar()]) -> erlynum:nvector().
%% @equiv from_list(List, [{dtype, auto}])
from_list(List) -> from_list(List, [{dtype, auto}]).

-spec from_list([erlynum:nscalar()], [erlynum:create_option()]) -> erlynum:nvector().
%% @doc Returns a vector from the given list items.
from_list(List, Options) -> ?WE(erlynum_nif:nvector_from_list(List, Options)).

-spec to_list(NVector :: erlynum:nvector()) ->  [ erlynum:nscalar() ].
%% @equiv to_list(NVector, noconvert)
to_list(NVector) -> to_list(NVector, noconvert).

-spec to_list(
    NVector         ::  erlynum:nvector(),
    Convert         ::  erlynum:convert_option()
) ->  [ erlynum:nscalar() ].
%% @doc Returns an Erlang list from the given vector value.
%% @param Convert specifies target data type.
to_list(NVector, Convert) -> ?WE(erlynum_nif:nvector_to_list(NVector, Convert)).

-spec get(
    NVector         ::  erlynum:nvector(),
    Index           ::  non_neg_integer()
) -> erlynum:nscalar().
%% @equiv get(NVector, Index, noconvert)
get(NVector, Index) -> get(NVector, Index, noconvert).

-spec get(
    NVector         ::  erlynum:nvector(),
    Index           ::  non_neg_integer(),
    Convert         ::  erlynum:convert_option()
) -> erlynum:nscalar().
%% @doc Returns a scalar element value from the `Vector'.
%% Note that `Index' starts from 0, but not 1.
%% It is possible to specify `Convert' option to cast scalar value into integer, real of complex value.
get(NVector, Index, Convert) -> ?WE(erlynum_nif:nvector_get(NVector, Index, Convert)).

-spec range(erlynum:nscalar()) -> erlynum:nvector().
%% @equiv range(0, Stop, 1, [])
range(Stop) -> range(0, Stop, 1, []).

-spec range(
    Start           ::  erlynum:nscalar(),
    Options_Stop    ::  [ erlynum:create_option() | erlynum:range_option() ]
                    |   erlynum:nscalar()
) -> erlynum:nvector().
%% @doc Returns a vector filled by values withing from `Start' and step `1'.
%%
%% If last argument is a {@link erlynum:scalar()} value, it specifies range `Stop'.
%% If last argument is a list, it specifies create options.
%%
%% @see range/4
range(Stop, Options) when is_list(Options) ->
    range(0, Stop, 1, Options);
range(Start, Stop) ->
    range(Start, Stop, 1, []).

-spec range(
    Start           ::  erlynum:nscalar(),
    Stop            ::  erlynum:nscalar(),
    Step_Options    ::  [ erlynum:create_option() | erlynum:range_option() ]
                    |   erlynum:nscalar()
) -> erlynum:nvector().
%% @doc Returns a vector filled by values withing the given range `[Start, Stop)'.
%%
%% If last argument is a {@link erlynum:scalar()} value, it specifies increment on next value.
%% If last argument is a list, it specifies create options.
%%
%% @see range/4
range(Start, Stop, Step_Options) when is_list(Step_Options) ->
    range(Start, Stop, 1, Step_Options);
range(Start, Stop, Step_Options) ->
    range(Start, Stop, Step_Options, []).

-spec range(
    Start           ::  erlynum:nscalar(),
    Stop            ::  erlynum:nscalar(),
    Step            ::  erlynum:nscalar(),
    Options         :: [ erlynum:create_option() | erlynum:range_option() ]
) -> erlynum:nvector().
%% @doc Returns a vector filled by values withing the given range `[Start, Stop)' and increment `Step'.
range(Start, Stop, Step, Options) -> ?WE(erlynum_nif:nvector_range(Start, Stop, Step, Options)).

-spec linspace(
    Start           ::  erlynum:nscalar(),
    Stop            ::  erlynum:nscalar(),
    Count           ::  non_neg_integer()
) -> erlynum:nvector().
%% @equiv linspace(Start, Stop, Count, [])
linspace(Start, Stop, Count) -> linspace(Start, Stop, Count, []).

-spec linspace(
    Start           ::  erlynum:nscalar(),
    Stop            ::  erlynum:nscalar(),
    Count           ::  non_neg_integer(),
    Options         ::  [ erlynum:create_option() | erlynum:range_option() ]
) -> erlynum:nvector().
%% @doc Returns a vector filled by the specified `Count' of values  withing the given range `[Start, Stop)'.
linspace(Start, Stop, Count, Options) -> ?WE(erlynum_nif:nvector_linspace(Start, Stop, Count, Options)).

-spec logspace(
    Start           ::  erlynum:nscalar(),
    Stop            ::  erlynum:nscalar(),
    Count           ::  non_neg_integer()
) -> erlynum:nvector().
%% @equiv logspace(Stat, Stop, Count, 10, [])
logspace(Start, Stop, Count) ->
    logspace(Start, Stop, Count, 10.0).

-spec logspace(
    Start           ::  erlynum:nscalar(),
    Stop            ::  erlynum:nscalar(),
    Count           ::  erlynum:nscalar(),
    Base_Options    ::  [ erlynum:create_option() | erlynum:range_option() ]
                    |   erlynum:nscalar()
) -> erlynum:nvector().
%% @doc Returns a vector filled by the specified `Count' of powers withing the given range `[Start, Stop)'.
%%
%% If last argument is a {@link erlynum:scalar()} value, it specifies `Base'.
%% If last argument is a list, it specifies create options.
%%
%% @see logspace/5
logspace(Start, Stop, Count, Base_Options) when is_number(Base_Options) or is_tuple(Base_Options) ->
    logspace(Start, Stop, Count, Base_Options, []);
logspace(Start, Stop, Count, Base_Options) when is_list(Base_Options) ->
    logspace(Start, Stop, Count, 10.0, Base_Options).

-spec logspace(
    Start           ::  erlynum:nscalar(),
    Stop            ::  erlynum:nscalar(),
    Count           ::  non_neg_integer(),
    Base            ::  erlynum:nscalar(),
    Options         ::  [ erlynum:create_option() | erlynum:range_option() ]
) -> erlynum:nvector().
%% @doc Returns a vector filled by the specified `Count' of powers of `Base' withing the given range `[Start, Stop)'.
logspace(Start, Stop, Count, Base, Options) -> ?WE(erlynum_nif:nvector_logspace(Start, Stop, Count, Base, Options)).

-spec geomspace(
    Start           ::  erlynum:nscalar(),
    Stop            ::  erlynum:nscalar(),
    Count           ::  non_neg_integer()
) -> erlynum:nvector().
%% @equiv geomspace(Start, Stop, Count, [])
geomspace(Start, Stop, Count) -> geomspace(Start, Stop, Count, []).

-spec geomspace(
    Start           ::  erlynum:nscalar(),
    Stop            ::  erlynum:nscalar(),
    Count           ::  non_neg_integer(),
    Options         ::  [ erlynum:create_option() | erlynum:range_option() ]
) -> erlynum:nvector().
%% @doc Returns a vector filled by specified `Count' of geometric progression values.
geomspace(Start, Stop, Count, Options) -> ?WE(erlynum_nif:nvector_geomspace(Start, Stop, Count, Options)).

-spec copy(X :: erlynum:nvector()) -> erlynum:nvector().
%% @equiv copy(X, [])
copy(X) -> copy(X, []).

-spec copy(X :: erlynum:nvector(), Options :: [erlynum:create_option()]) -> erlynum:nvector().
%% @doc Returns a packed copy of vector view `Y ← αX'.
%% This allow Erlang garbage collector to release referenced binary.
%% The function is useless in case of vector view source (some big matrix, for example)
%% stored in the same process stack, due to referenced source binary still alive and will not be collected.
%%
%% Generic cases to use this function are:
%%   * to pass return value to another process or node if source vector generated from another big array-like object;
%%   * to convert internal data type while copying.
copy(X, Options) -> ?WE(erlynum_nif:nvector_copy(X, Options)).

-spec scale(
    Y               ::  erlynum:nvector(),
    Alpha           ::  erlynum:nscalar()
) -> erlynum:nvector().
%% @equiv scale(Y, Alpha, [])
scale(Y, Alpha) -> scale(Y, Alpha, []).

-spec scale(
    Y               ::  erlynum:nvector(),
    Alpha           ::  erlynum:nscalar(),
    Options         ::  [ erlynum:create_option() ]
) -> erlynum:nvector().
%% @doc Returns a scaled copy `Y ← αY' of the vector of the same size and the specified data type.
scale(Y, Alpha, Options) -> ?WE(erlynum_nif:nvector_scale(Y, Alpha, Options)).

-spec sum(Y :: erlynum:nvector(), X :: erlynum:nvector()) -> erlynum:nvector().
%% @equiv sum(Y, X, [])
sum(NVectorY, NVectorX) -> axpy(NVectorX, NVectorY, 1).

-spec sum(
    Y               ::  erlynum:nvector(),
    X               ::  erlynum:nvector(),
    Options         ::  [ erlynum:create_option() ]
) -> erlynum:nvector().
%% @doc Returns an element-wise sum of two vectors `Y ← X + Y'.
%%
%% If vector sizes differs, the result is a vector of greater size, leaving
%% rest elements unchanged.
%%
%% If data types differs, the `Y' vector (first argument) data type is used until options provided.
sum(Y, X, Options) -> axpy(X, Y, 1, Options).

-spec axpy(
    Y               ::  erlynum:nvector(),
    X               ::  erlynum:nvector(),
    Alpha           ::  erlynum:nscalar()
) -> erlynum:nvector().
%% @equiv axpy(Y, X, Alpha, [])
axpy(Y, X, Alpha) -> axpy(Y, X, Alpha, []).

-spec axpy(
    Y               ::  erlynum:nvector(),
    X               ::  erlynum:nvector(),
    Alpha           ::  erlynum:nscalar(),
    Options         ::  [ erlynum:create_option() ]
) -> erlynum:nvector().
%% @doc Returns a vector, each element calculated in form `Y ← αX + Y'
axpy(Y, X, Alpha, Options) ->
    ?WE(erlynum_nif:nvector_axpy(Y, X, Alpha, Options)).

-spec dot(
    Y               :: erlynum:nvector(),
    X               :: erlynum:nvector()
) -> erlynum:nscalar().
%% @equiv dot(Y, X, [])
dot(Y, X) -> dot(Y, X, []).

-spec dot(
    Y               :: erlynum:nvector(),
    X               :: erlynum:nvector(),
    Options         :: [ erlynum:create_option() | {conjuated, boolean()}]
) -> erlynum:nscalar().
%% @doc Returns the inner product of two vectors.
dot(Y, X, Options) -> ?WE(erlynum_nif:nvector_dot(Y, X, Options)).

-spec asum(X :: erlynum:nvector()) -> erlynum:nscalar().
%% @equiv asum(X, [])
asum(X) -> asum(X, []).

-spec asum(
    X               :: erlynum:nvector(),
    Options         :: [ erlynum:create_option() ]
) -> erlynum:nscalar().
%% @doc Returns the sum of magnitudes of elements of a real vector, or the sum of magnitudes
%% of the real and imaginary parts of elements of a complex vector.
asum(X, Options) -> ?WE(erlynum_nif:nvector_asum(X, Options)).

-spec iamax(X :: erlynum:nvector()) -> non_neg_integer() | undefined.
%% @doc Returns the lowest index of vector element that has the largest absolute value.
iamax(X) -> ?WE(erlynum_nif:nvector_iamax_iamin(X, iamax)).

-spec iamin(X :: erlynum:nvector()) -> non_neg_integer() | undefined.
%% @doc Returns the lowest index of vector element that has the lowest absolute value.
iamin(X) -> ?WE(erlynum_nif:nvector_iamax_iamin(X, iamin)).

-spec nrm2(X :: erlynum:nvector()) -> erlynum:nscalar().
%% @equiv nrm2(X, [])
nrm2(X) -> nrm2(X, []).

-spec nrm2(
    X               :: erlynum:nvector(),
    Options         :: [ erlynum:create_option() ]
) -> erlynum:nscalar().
%% @doc Returns the Euclidean norm of vector.
%%
%% The function is named `nrm2' to match corresponding function name in BLAST.
%% It is better to use `euclidean_norm' function for better readability.
nrm2(X, Options) -> ?WE(erlynum_nif:nvector_nrm2(X, Options)).

-spec euclidean_norm(X :: erlynum:nvector()) -> erlynum:nscalar().
%% @doc equiv nrm2(X, [])
euclidean_norm(X) -> nrm2(X).

-spec euclidean_norm(
    X               :: erlynum:nvector(),
    Options         :: [ erlynum:create_option() ]
) -> erlynum:nscalar().
%% @equiv nrm2(X, Options)
euclidean_norm(X, Options) -> nrm2(X, Options).



% --------------- EUnit testing functions

zeros_empty_test() ->
    NVector = zeros(0),
    {RowsCount, ColsCount} = NVector#nvector.shape,
    TotalCount = RowsCount * ColsCount,
    ?assertEqual(0, TotalCount),
    BinData = NVector#nvector.data,
    ?assertEqual(0, byte_size(BinData)),
    ?assertEqual(0, NVector#nvector.view#view_params.size).

zeros_nonempty_test() ->
    Zeros = zeros(5),
    ?assertEqual([0.0, 0.0, 0.0, 0.0, 0.0], to_list(Zeros)),
    ?assertEqual([0, 0, 0, 0, 0], to_list(Zeros, integer)),
    ComplexZeros = zeros(5, [{dtype, c}]),
    ?assertEqual([{0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}], to_list(ComplexZeros)),
    ?assertEqual([0, 0, 0, 0, 0], to_list(ComplexZeros, integer)).

ones_empty_test() ->
    NVector = ones(0),
    {RowsCount, ColsCount} = NVector#nvector.shape,
    TotalCount = RowsCount * ColsCount,
    ?assertEqual(0, TotalCount),
    BinData = NVector#nvector.data,
    ?assertEqual(0, byte_size(BinData)),
    ?assertEqual(0, NVector#nvector.view#view_params.size).

ones_nonempty_test() ->
    Ones = ones(5),
    ?assertEqual([1.0, 1.0, 1.0, 1.0, 1.0], to_list(Ones)),
    ?assertEqual([1, 1, 1, 1, 1], to_list(Ones, integer)),
    ComplexOnes = ones(5, [{dtype, c}]),
    ?assertEqual([{1.0,0.0}, {1.0,0.0}, {1.0,0.0}, {1.0,0.0}, {1.0,0.0}], to_list(ComplexOnes)),
    ?assertEqual([1, 1, 1, 1, 1], to_list(ComplexOnes, integer)).

full_empty_test() ->
    NVector = full(0, 123),
    {RowsCount, ColsCount} = NVector#nvector.shape,
    TotalCount = RowsCount * ColsCount,
    ?assertEqual(0, TotalCount),
    BinData = NVector#nvector.data,
    ?assertEqual(0, byte_size(BinData)),
    ?assertEqual(0, NVector#nvector.view#view_params.size).

full_nonempty_test() ->
    Full = full(5, 123),
    ?assertEqual([123.0, 123.0, 123.0, 123.0, 123.0], to_list(Full)),
    ?assertEqual([123, 123, 123, 123, 123], to_list(Full, integer)),
    ComplexFull = full(3, {123,456}),
    ?assertEqual([{123.0,456.0}, {123.0,456.0}, {123.0,456.0}], to_list(ComplexFull)),
    ?assertEqual([123.0, 123.0, 123.0], to_list(ComplexFull, real)),
    ?assertEqual([123, 123, 123], to_list(ComplexFull, integer)).

asum_real_test() ->
    Zeros = zeros(5),
    ?assertEqual(0.0, asum(Zeros)),
    Ones = ones(5),
    ?assertEqual(5.0, asum(Ones)),
    PositiveRange = from_list([1,2,3,4]),
    VariativeRange = from_list([-1,2,-3,4]),
    ?assertEqual(1.0+2.0+3.0+4.0, asum(PositiveRange)),
    ?assertEqual(1.0+2.0+3.0+4.0, asum(VariativeRange)).

asum_complex_test() ->
    Zeros = zeros(5, [{dtype, z}]),
    ?assertEqual(0.0, asum(Zeros)),
    Ones = ones(5, [{dtype, z}]),
    ?assertEqual(5.0, asum(Ones)),
    PositiveRange = from_list([{1.0,2.0}, {3.0,4.0}]),
    VariativeRange = from_list([{-1.0,2.0}, {-3.0,4.0}]),
    ?assertEqual(1.0+2.0+3.0+4.0, asum(PositiveRange)),
    ?assertEqual(1.0+2.0+3.0+4.0, asum(VariativeRange)).

iamax_test() ->
    ?assertEqual(undefined, iamax(from_list([]))),
    %                                0  1  2  3  4  5  6  7
    ?assertEqual(4, iamax(from_list([1,-2, 3, 2,-5, 1, 4, 5]))).

iamin_test() ->
    ?assertEqual(undefined, iamin(from_list([]))),
    %                                  0  1  2   3  4  5  6
    ?assertEqual(4, iamin(from_list([ -2, 3, 2, -5, 1, 4, 5]))).

euclidean_norm_test() ->
    ?assertEqual(0.0, euclidean_norm(zeros(5))).