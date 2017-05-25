%% @author          Victor Yacovlev <v.yacovlev@gmail.com>
%% @copyright       2017 Victor Yacovlev
%% @doc Type definitions for entire library to be used by Dialyzer.
-module(erlynum).
-include("erlynum.hrl").
-export_type([
    precision/0, dtype/0, shape/0,
    complex/0, nvector/0, nmatrix/0, nscalar/0, narray/0,
    driver_backend/0, create_option/0, range_option/0, convert_option/0, view_params/0
]).


%% @type precision() = single | double.
%% Preferred floating point precision.
-type precision()       ::  single | double.

%% @type dtype() = s | d | c | z.
%% Internal data type for one element.
%% Atom `s' is single precision real value, `d' is double precision real,
%% `c' is single precision complex value, `z' is double precision complex.
-type dtype()           ::  s | d | c | z.

%% @type shape() = non_neg_integer() | { Rows :: non_neg_integer(), Columns :: non_neg_integer()}.
%% Array shape. Just a number of items for vector, or `{Rows, Columns}' pair for matrix.
-type shape()           ::  non_neg_integer() | {non_neg_integer(), non_neg_integer()}.

%% @type complex() = { Real :: number(), Imaginary :: number()}.
%% Complex value representation.
%% The first item is a `Real' part of complex number, the second is a `Imaginary'.
-type complex()         ::  {number(), number()}.

% -record(view_params, {
%     offset              ::  non_neg_integer(),
%     size                ::  non_neg_integer(),
%     increment           ::  integer()
% }).

%% @type view_params() = #view_params{}.
%% Internal view on the array data for the each dimension.
%% Field `size' is a number of items in view, `offset' is an index of first view element in the array,
%% field `increment' is a distance between elements (might be negative for reverse or transposed views).
-type view_params()     ::  #view_params{}.

%% @type nvector() = #nvector{}.
%% One-dimensional vector.
-type nvector()         ::  #nvector{}.

%% @type nmatrix() = #nmatrix{}.
%% Two-dimensional matrix.
-type nmatrix()         ::  #nmatrix{}.

-type driver_backend()  ::  blas | atlas | mkl.

-type create_option()   ::  {dtype, auto | dtype()}
                        |   {precision, precision()}.

-type range_option()    ::  {endpoint, boolean()}.

-type convert_option()  ::  noconvert | integer | real | complex.

-type narray()          ::  nvector() | nmatrix().

-type nscalar()         ::  number() | complex().