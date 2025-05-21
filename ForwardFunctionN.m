function [t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan)

%Compute forward solution to initial value problem to get tumor volume
%curves for each sensitivity value, then take their weighted sum

%INPUTS: 
    %sgrid: mesh of sensitivity distribution
    %sprobs: sensitivity distribution
    %rho: growth rate
    %k: death rate
    %y0: initial scaled cell count
    %tspan: time vector

%OUTPUTS: 
    %t: time vector returned by ode45
    %cmat: matrix of solution curves for each sensitivity value
    %weightedsol: weighted sum of solution curves, representing aggregate
    %tumor volume 

    cmat=zeros(length(tspan),length(sgrid));
    for i=1:length(sgrid)
        %solve forward problem for each s value
        [t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) -k*sgrid(i)*csol, tspan, y0);
        cmat(:,i)=csol;
    end
    
    %compute weighted sum
    weightedsol=(sprobs./(sum(sprobs)))*cmat';
 end 