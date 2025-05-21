function[t,c]=ForwardFnFindRho(rho,tspan,IC)

%Solution of the forward problem used in finding the minimum and maximum 
%growth rate rho. 

%INPUTS:
    %rho: growth rate 
    %tspan: time vector from data
    %IC: initial condition from data

%OUTPUTS: 
    %t: time vector given by ode45
    %c: aggregated tumor volume given by ode45

%solve forward problem with only growth term 
[t,c]=ode45(@(t,c) rho*c*(1-c),tspan,IC );
