function [error_call,error_put] = pricing_error(call_BSM,call,put_BSM,put)

% ERROR RELATIVE TO THE BLACK-SCHOLES-MERTON FUNCTION

error_call = abs((call - call_BSM)/call_BSM);
error_put = abs((put - put_BSM)/put_BSM);

end

