%find the weighted sum of squared differences between data and model
function[err,bestapprox]=PEerrorfn(s, gencmatC, data,weights)

bestapprox=(gencmatC*s)'; 

err=sum((weights.^(1)).*((bestapprox - data)).^2);

end