function[t,c]=ForwardFnFindK(k,tspan,IC)

[t,c]=ode45(@(t,c) -k*c,tspan,IC ); 
