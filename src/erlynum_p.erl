%% @author          Victor Yacovlev <v.yacovlev@gmail.com>
%% @copyright       2017 Victor Yacovlev
%% @private
-module(erlynum_p).

-export([wrap_error/1]).

wrap_error({error, Error}) -> erlang:error(Error);
wrap_error(NifRetValue) -> NifRetValue.
