%Function that takes in weightedsol and constdata, as well as r, and finds 
% the sum of squared differences

function[err,bestapprox]=PEerrorfn(s, gencmatC, data,weights)

bestapprox=(gencmatC*s)'; %bestapprox gets weightedsol by multiplying s*csol

err=sum((weights.^(1)).*((bestapprox - data)).^2);

end