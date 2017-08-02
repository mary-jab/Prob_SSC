function [id, minError, error]= CRC_IDcheck(D,class_pinv_M,y,Dlabels, maxCluster) % maryam add error
%------------------------------------------------------------------------
% CRC_RLS classification function
if nargin <5
    maxCluster = 0;
end
maxCluster = max(maxCluster, max(Dlabels));
coef         =  class_pinv_M*y;
for ci = 1:maxCluster
    coef_c   =  coef(Dlabels==ci);
    Dc       =  D(:,Dlabels==ci);
    error(ci) = norm(y-Dc*coef_c,2)^2/sum(coef_c.*coef_c);
end

minError = min(error);
index      =  find(error==minError);
id         =  index(1);
end