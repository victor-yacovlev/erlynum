-ifndef(ERLYNUM_HRL).
-define(ERLYNUM_HRL, 1).

-record(view_params, {
    offset              ::  non_neg_integer(),
    size                ::  non_neg_integer(),
    increment           ::  integer()
}).

-record(nvector, {
    view                ::  erlynum:view_params(),
    shape               ::  erlynum:shape(),
    dtype       = d     ::  erlynum:dtype(),
    data        = <<>>  ::  binary()
}).

-record(nmatrix, {
    view                ::  { erlynum:view_params(), erlynum:view_params() },
    shape               ::  erlynum:shape(),
    dtype       = d     ::  erlynum:dtype(),
    data        = <<>>  ::  binary()
}).

-endif.