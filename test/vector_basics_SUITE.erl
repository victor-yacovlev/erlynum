-module(vector_basics_SUITE).

-include_lib("eunit/include/eunit.hrl").
-include("erlynum.hrl").

zeros_empty_test() ->
    NVector = nvector:zeros(0),
    {RowsCount, ColsCount} = NVector#nvector.shape,
    TotalCount = RowsCount * ColsCount,
    ?assertEqual(0, TotalCount),
    BinData = NVector#nvector.data,
    ?assertEqual(0, byte_size(BinData)),
    ?assertEqual(0, NVector#nvector.view#view_params.size).

zeros_nonempty_test() ->
    Zeros = nvector:zeros(5),
    ?assertEqual([0.0, 0.0, 0.0, 0.0, 0.0], nvector:to_list(Zeros)),
    ?assertEqual([0, 0, 0, 0, 0], nvector:to_list(Zeros, integer)),
    ComplexZeros = nvector:zeros(5, [{dtype, c}]),
    ?assertEqual([{0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}], nvector:to_list(ComplexZeros)),
    ?assertEqual([0, 0, 0, 0, 0], nvector:to_list(ComplexZeros, integer)).

ones_empty_test() ->
    NVector = nvector:ones(0),
    {RowsCount, ColsCount} = NVector#nvector.shape,
    TotalCount = RowsCount * ColsCount,
    ?assertEqual(0, TotalCount),
    BinData = NVector#nvector.data,
    ?assertEqual(0, byte_size(BinData)),
    ?assertEqual(0, NVector#nvector.view#view_params.size).

ones_nonempty_test() ->
    Ones = nvector:ones(5),
    ?assertEqual([1.0, 1.0, 1.0, 1.0, 1.0], nvector:to_list(Ones)),
    ?assertEqual([1, 1, 1, 1, 1], nvector:to_list(Ones, integer)),
    ComplexOnes = nvector:ones(5, [{dtype, c}]),
    ?assertEqual([{1.0,0.0}, {1.0,0.0}, {1.0,0.0}, {1.0,0.0}, {1.0,0.0}], nvector:to_list(ComplexOnes)),
    ?assertEqual([1, 1, 1, 1, 1], nvector:to_list(ComplexOnes, integer)).

full_empty_test() ->
    NVector = nvector:full(0, 123),
    {RowsCount, ColsCount} = NVector#nvector.shape,
    TotalCount = RowsCount * ColsCount,
    ?assertEqual(0, TotalCount),
    BinData = NVector#nvector.data,
    ?assertEqual(0, byte_size(BinData)),
    ?assertEqual(0, NVector#nvector.view#view_params.size).

full_nonempty_test() ->
    Full = nvector:full(5, 123),
    ?assertEqual([123.0, 123.0, 123.0, 123.0, 123.0], nvector:to_list(Full)),
    ?assertEqual([123, 123, 123, 123, 123], nvector:to_list(Full, integer)),
    ComplexFull = nvector:full(3, {123,456}),
    ?assertEqual([{123.0,456.0}, {123.0,456.0}, {123.0,456.0}], nvector:to_list(ComplexFull)),
    ?assertEqual([123.0, 123.0, 123.0], nvector:to_list(ComplexFull, real)),
    ?assertEqual([123, 123, 123], nvector:to_list(ComplexFull, integer)).