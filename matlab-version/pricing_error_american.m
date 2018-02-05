function [error_call_American,error_put_American] = pricing_error_american(call_BIN_American,call_American,put_BIN_American,put_American)

% ERROR RELATIVE TO THE AMERICAN BINOMIAL FUNCTION

error_call_American = abs((call_American - call_BIN_American)/call_BIN_American);
error_put_American = abs((put_American - put_BIN_American)/put_BIN_American);

end