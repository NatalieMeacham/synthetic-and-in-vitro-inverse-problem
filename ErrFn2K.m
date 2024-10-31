function[err]=ErrFn2K(k, data,tspan,IC,weights)

rhoval=0.0682;

%bestapprox=(gencmatC*s)'; %bestapprox gets weightedsol by multiplying s*csol
%[~,c]=ode45(@(t,c) -k*c,tspan,IC );
%[~,c]=ode45(@(t,c) -k*c*rhoval,tspan,IC );
[~,c]=ode45(@(t,c) -k*c*(1-c),tspan,IC );

err=sum((weights).*((c' - data)).^2);
    %sum((weights).*((bestapprox - data)).^2);