
% Log-likelihoood for Probit function

function val = ProbitLogLike(b,y,x)
val1 = y'*log(normcdf(x*b));
val0 = (1-y)'*log(1-normcdf(x*b));
val = -1.*(sum(val1 + val0));

end
