
function[err,bestapprox]=PEerrorfn(s, gencmatC, data,weights)

%find the weighted sum of squared differences between data and model

%INPUTS:
    %s: vector of sensitivity values
    %genmatc: matrix of forward solutions for each s value
    %data: synthetic or in vitro data vector
    %weights: weights from gls 

%OUTPUTS: 
    %err: value of summed error
    %bestapprox: approximation of aggregate tumor volume using genmatC and
    %s

bestapprox=(gencmatC*s)'; 

err=sum((weights.^(1)).*((bestapprox - data)).^2);

end