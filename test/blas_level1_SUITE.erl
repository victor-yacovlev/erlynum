-module(blas_level1_SUITE).

-include_lib("eunit/include/eunit.hrl").

asum_real_test() ->
    Zeros = nvector:zeros(5),
    ?assertEqual(0.0, nvector:asum(Zeros)),
    Ones = nvector:ones(5),
    ?assertEqual(5.0, nvector:asum(Ones)),
    PositiveRange = nvector:from_list([1,2,3,4]),
    VariativeRange = nvector:from_list([-1,2,-3,4]),
    ?assertEqual(1.0+2.0+3.0+4.0, nvector:asum(PositiveRange)),
    ?assertEqual(1.0+2.0+3.0+4.0, nvector:asum(VariativeRange)).

asum_complex_test() ->
    Zeros = nvector:zeros(5, [{dtype, z}]),
    ?assertEqual(0.0, nvector:asum(Zeros)),
    Ones = nvector:ones(5, [{dtype, z}]),
    ?assertEqual(5.0, nvector:asum(Ones)),
    PositiveRange = nvector:from_list([{1.0,2.0}, {3.0,4.0}]),
    VariativeRange = nvector:from_list([{-1.0,2.0}, {-3.0,4.0}]),
    ?assertEqual(1.0+2.0+3.0+4.0, nvector:asum(PositiveRange)),
    ?assertEqual(1.0+2.0+3.0+4.0, nvector:asum(VariativeRange)).

iamax_test() ->
    ?assertEqual(undefined, nvector:iamax(nvector:from_list([]))),
    %                                                0  1  2  3  4  5  6  7
    ?assertEqual(4, nvector:iamax(nvector:from_list([1,-2, 3, 2,-5, 1, 4, 5]))).

iamin_test() ->
    ?assertEqual(undefined, nvector:iamin(nvector:from_list([]))),
    %                                                  0  1  2   3  4  5  6
    ?assertEqual(4, nvector:iamin(nvector:from_list([ -2, 3, 2, -5, 1, 4, 5]))).

euclidean_norm_test() ->
    ?assertEqual(0.0, nvector:euclidean_norm(nvector:zeros(5))).