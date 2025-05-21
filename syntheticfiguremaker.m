%Script to run synthetic data inverse problem for different initial
%distributions, noise levels, and time meshes.
%Generates figures 2-4 in the main paper as well as figures 1-6 in the supplementary material.  

disttypes = {'Normal','Bigaussian','Uniform','OnePoint','TwoPoints'};
noiselevels = [0, 0.01, 0.02, 0.05];
timepoints = [5, 10, 25, 100];

%% Run IP for given distribution, noise level, and time mesh
%To generate supplementary figure 1, use the Bigaussian distribution option
%and manually change the definition of that distribution in DistFn2

[sgrid,sprobs, rsgrid,optweightfromAIC,t,propdata,sweightedsol,cdfS,cdfR]= GLSInverseScriptNPaper('TwoPoints', 100,0.0);

%% Noiseless synthetic data for all five distributions (figure 2 in paper)
%Also gives AIC curves from supplementary figure 2 and non-scaled recovered
%results from supplementary figure 3

for kk=1:5 
    GLSInverseScriptNPaper(disttypes{kk},timepoints(end),noiselevels(1))
end

%% Comparing time points (figure 3 in main paper)
%Also includes all other distributions (supplementary figure 4)

for kk =1:length(disttypes)
    [sgrid,sprobs,rsgrid5,optweight5, t, propdata5, weightedsol5,cdfS, cdfR5] = GLSInverseScriptNPaper(disttypes{kk},timepoints(1),noiselevels(1));
    [~,~,rsgrid10,optweight10, t10, propdata10, weightedsol10,~, cdfR10] = GLSInverseScriptNPaper(disttypes{kk},timepoints(2),noiselevels(1));
    [~,~,rsgrid25,optweight25, t50, propdata25, weightedsol25,~, cdfR25] = GLSInverseScriptNPaper(disttypes{kk},timepoints(3),noiselevels(1));
    [~,~,rsgrid100,optweight100,t100, propdata100, weightedsol100,~, cdfR100] = GLSInverseScriptNPaper(disttypes{kk},timepoints(4),noiselevels(1));

    % create figures combining all time meshes
    timestepplots(sprobs,optweight100,optweight25,optweight10,optweight5,propdata100, propdata25,propdata10,propdata5,weightedsol100,weightedsol25,weightedsol10,weightedsol5,disttypes{kk},t100,propdata100, noiselevels(1))
end

%% Noisy data recovery (figure 4 in paper)
%Change noise level to get figures 5 and 6 from supplementary material

for kk=1:length(disttypes) 
    GLSInverseScriptNPaper(disttypes{kk},100,noiselevels(3))
end