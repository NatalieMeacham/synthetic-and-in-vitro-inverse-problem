%Function that takes in weightedsol and constdata, as well as r, and finds 
% the sum of squared differences

%note that for gls, will need to add weights to this function

%function[err,bestapprox]=PEerrorfn(s, gencmat, constdata,weights)
%function[err,bestapprox]=PEerrorfn(s, gencmatC, constdata,weights)
function[err,bestapprox]=PEerrorfn(s, gencmatC, data,weights)

bestapprox=(gencmatC*s)'; %bestapprox gets weightedsol by multiplying s*csol

%err=sum(weights.*((bestapprox - constdata).^2));
err=sum((weights.^(1)).*((bestapprox - data)).^2);
%err=sum((weights.^(-2)).*((bestapprox - data)).^2);
%err=sum((weights.^(-2)).*((bestapprox - data)').^2)

%note cmat shape is t points # of rows and s points # of cols. so 

%it seems like r does not need to be here if we do not need to reshape?