function [sprobs] = DistFn2(disttype,sgrid,a,b)

%FUTURE UPDATE: do "BigaussianLow" and "BigaussianHigh" for two
%different bigaussian dists

if string(disttype) == 'Normal' %this should be cts
    mu = .5; 
    sigma = .09;
    sprobs = pdf('Normal',sgrid,mu,sigma); 
    %sprobs = sprobs/(sum(sprobs));
% elseif string(disttype) == 'CtsNormal'
%     mu = .5; 
%     sigma = .09;
%     %sprobs = pdf('Normal',sgrid,mu,sigma); 
%     %normal = @(x) (1/(sigma*sqrt(2*pi)))*exp(-.5*((x - mu)/sigma).^2);
%     %sprobs=normal(sgrid)
%     sprobs = normpdf(sgrid,mu,sigma); 
 elseif string(disttype) == 'Bigaussian' %this should be cts 
     %mu1=.3;
     %mu2=.7; 
     mu1 = 0.15;
     mu2 = 0.85;
     % sigma1=.05;
     % sigma2=.1;
     % weight1=.1;
     % weight2=.9;
     sigma1=.1; %original 
     sigma2=.05; %original
     % sigma1 = 0.1
     % sigma2 = 0.1 %evening version
     weight1=.9; %typically use this pair
     weight2=.1;
     % weight1=.1;
     % weight2=.9;
      weight1 = 0.5; %paper new version
       weight2 = 0.5;
     %  %bigaussian = @(x) (1/(sigma1*sqrt(2*pi)))*exp(-.5*((x - mu1)/sigma1).^2) + (1/(sigma2*sqrt(2*pi)))*exp(-.5*((x - mu2)/sigma2).^2);
     bigaussian = @(x) weight1*(1/(sigma1*sqrt(2*pi)))*exp(-.5*((x - mu1)/sigma1).^2) + weight2*(1/(sigma2*sqrt(2*pi)))*exp(-.5*((x - mu2)/sigma2).^2);
     sprobs=bigaussian(sgrid);
     %sprobs = sprobs/(sum(sprobs));
elseif string(disttype) == 'Uniform'
     sprobs = pdf('Uniform',sgrid, a, b);
     %sprobs = sprobs/(sum(sprobs));
elseif string(disttype) == 'OnePoint'
    pointval=0.4;
    %point=pointval*(length(sgrid) - 1);
    %point=(sgrid(pointval*100 + 1));
    sprobs=zeros(length(sgrid),1);
    sprobs(pointval*100)=1;
    sprobs=sprobs';
    sprobs = sprobs/(sum(sprobs));
elseif string(disttype) == 'TwoPoints'
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
    disp('Choose a distribution that works.')
end

%Normalize sprobs to sum to 1 for whichever distribution:

end