function[err]=ErrorFn2(k, data,tspan,IC,weights)

%bestapprox=(gencmatC*s)'; %bestapprox gets weightedsol by multiplying s*csol
[~,c]=ode45(@(t,c) -k*c,tspan,IC );
%[c]=ForwardFunctionFindRho(rho, data,tspan,IC)

err=sum((weights).*((c' - data)).^2);
    %sum((weights).*((bestapprox - data)).^2);