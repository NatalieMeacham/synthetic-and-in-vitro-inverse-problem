function [meanmatnorm]=ReplaceNaNsFn(concvec,meanmatnorm)
%test comment 

%replace any NaNs row value with 1s
TFvec=zeros(length(concvec),1);
for a=concvec
    TFvec(a)=anynan(meanmatnorm(a,:))
    %TFvec=[TFvec, TF]
end

for a=concvec
    if TFvec(a)==1
        meanmatnorm(a,:)=1;
        disp('updated meanmatnorm with ones')
    else 
        disp('chill')
    end
end