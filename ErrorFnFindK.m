function[err]=ErrorFnFindK(k, data,tspan,IC,weights)

%Error function used in side inverse problem finding maximum death rate k

    %INPUTS: 
    %k: death rate being optimized
    %data: normalized cell count data (averaged specific concentration from
    %specific dataset)
    %tspan: time vector from dataset
    %IC: initial condition from data 
    %weights: GLS weights

    %OUTPUTS:
    %err: error between data and model 

%solve forward problem    
[~,c]=ode45(@(t,c) -k*c,tspan,IC );

%compute error
err=sum((weights).*((c' - data)).^2);