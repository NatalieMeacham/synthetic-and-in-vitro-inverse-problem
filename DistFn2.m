function [sprobs] = DistFn2(disttype,sgrid,a,b)

%Take in the sensitivity mesh and return a continuous or discrete 
%probability distribution determined by disttype.

%INPUTS
    %disttype: disttype: distribution type for initial sensitivity
    %distribution
    %sgrid: sensitivity mesh, determined by points, a, and b
    %a: lowest point of sensitivity mesh, always 0
    %b: highest point of sensitivity mesh, always 1

%OUTPUTS:
    %sprobs: vector indiciating p(s) for each s value in sgrid

if string(disttype) == 'Normal' %continuous Gaussian distribution
    mu = .5; 
    sigma = .09;
    sprobs = pdf('Normal',sgrid,mu,sigma); 
    sprobs = sprobs/(trapz(sgrid,sprobs));

 elseif string(disttype) == 'Bigaussian' %continuous bi-Gaussian distribution
     %main text figure 2 second row
     mu1=.3;
     mu2=.7; 
     sigma1=.05;
     sigma2=.1;
     weight1=.5; 
     weight2=.5;

     %appendix figure 10 first row (right bump 10% weight)
     % mu1=.3;
     % mu2=.7; 
     % sigma2=.05;
     % sigma1=.1;
     % weight2=.1; 
     % weight1=.9;

     %appendix figure 10 second row (left bump 10% weight)
     % mu1=.3;
     % mu2=.7; 
     % sigma2=.1;
     % sigma1=.05;
     % weight2=.9; 
     % weight1=.1;

     %appendix figure 10 third row (equally weighted bumps)
     % mu1=.3;
     % mu2=.7; 
     % sigma1=.05;
     % sigma2=.05;
     % weight1=.5; 
     % weight2=.5;

     %appendix figure 10 fourth row (equally weighted narrow bumps)
     % mu1=.3;
     % mu2=.7; 
     % sigma1=.01;
     % sigma2=.01;
     % weight1=.5; 
     % weight2=.5;

     %appendix figure 10 fifth row (equally weighted bumps far apart)
     % mu1=.15;
     % mu2=.85; 
     % sigma1=.05;
     % sigma2=.05;
     % weight1=.5; 
     % weight2=.5;

     bigaussian = @(x) weight1*(1/(sigma1*sqrt(2*pi)))*exp(-.5*((x - mu1)/sigma1).^2) + weight2*(1/(sigma2*sqrt(2*pi)))*exp(-.5*((x - mu2)/sigma2).^2);
     sprobs=bigaussian(sgrid);
     sprobs = sprobs/(trapz(sgrid,sprobs));

elseif string(disttype) == 'Uniform' %continuous uniform distribution
     sprobs = pdf('Uniform',sgrid, a, b);
     sprobs = sprobs/(trapz(sgrid,sprobs));

elseif string(disttype) == 'OnePoint' %discrete one-point distribution
    pointval=0.4; 
    sprobs=zeros(length(sgrid),1);
    sprobs(pointval*100)=1;
    sprobs=sprobs';
    sprobs = sprobs/(sum(sprobs));

elseif string(disttype) == 'TwoPoints' %discrete two-point distribution
    pointval1=0.3; 
    pointval2=0.8;
    prop1=0.4; 
    prop2=1 - prop1;
    sprobs=zeros(length(sgrid),1);
    sprobs(pointval1*100)=prop1;
    sprobs(pointval2*100)=prop2;
    sprobs=sprobs';
    sprobs = sprobs/(sum(sprobs));

else 
    disp('Choose an available distribution.')
end

end
