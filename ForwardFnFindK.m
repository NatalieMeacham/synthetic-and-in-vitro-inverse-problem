function[t,c]=ForwardFnFindK(k,tspan,IC)

%Solution of the forward problem used in finding the maximum death rate k. 

%INPUTS:
    %k: death rate 
    %tspan: time vector from data
    %IC: initial condition from data

%OUTPUTS: 
    %t: time vector given by ode45
    %c: aggregated tumor volume given by ode45

%solve forward problem with only death term
[t,c]=ode45(@(t,c) -k*c,tspan,IC ); 
