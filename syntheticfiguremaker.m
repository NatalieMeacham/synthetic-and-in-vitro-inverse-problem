%Synthetic Script

disttypes = {'Normal','Bigaussian','Uniform','OnePoint','TwoPoints'};
noiselevels = [0, 0.01, 0.02, 0.05];
timepoints = [5, 10, 25, 100];

%% Noiseless Synthetic Data (Fig 2 in main paper)

for kk=1:length(disttypes) 
    close all
    GLSInverseScriptNPaper(disttypes{kk},timepoints(end),noiselevels(1))
end

%% Comparing Time Points (Fig 3 in main paper)
for kk =1:length(disttypes)
    close all
    [sgrid,sprobs,rsgrid5,optweight5, t, propdata, weightedsol5,cdfS, cdfR5] = GLSInverseScriptNPaper(disttypes{kk},timepoints(1),noiselevels(1));
    [~,~,rsgrid10,optweight10, t10, propdata10, weightedsol10,~, cdfR10] = GLSInverseScriptNPaper(disttypes{kk},timepoints(2),noiselevels(1));
    [~,~,rsgrid25,optweight25, t50, propdata25, weightedsol25,~, cdfR25] = GLSInverseScriptNPaper(disttypes{kk},timepoints(3),noiselevels(1));
    [~,~,rsgrid100,optweight100,t100, propdata100, weightedsol100,~, cdfR100] = GLSInverseScriptNPaper(disttypes{kk},timepoints(4),noiselevels(1));

    % create figures
    timestepplots(sprobs,optweight100,optweight25,optweight10,propdata100, propdata25,propdata10,weightedsol100,weightedsol25,weightedsol10,disttypes{kk},t100,propdata100, noiselevels(1))
end
%% Noisy Cases              (Fig 4in main paper, 2 figs in Supp)
for kk=1:length(disttypes) 
    for jj =2:length(noiselevels)
        close all
        GLSInverseScriptNPaper(disttypes{kk},100,noiselevels(jj))
    end
end

