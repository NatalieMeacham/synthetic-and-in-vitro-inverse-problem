function[t,c]=ForwardFnFindRho(rho,tspan,IC)
[t,c]=ode45(@(t,c) rho*c*(1-c),tspan,IC );
