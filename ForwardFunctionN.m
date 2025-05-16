function [t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan)
    
    cmat=zeros(length(tspan),length(sgrid));
    for i=1:length(sgrid)
        [t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) -k*sgrid(i)*csol, tspan, y0);
        cmat(:,i)=csol;
    end
    
    weightedsol=(sprobs./(sum(sprobs)))*cmat';
 end 