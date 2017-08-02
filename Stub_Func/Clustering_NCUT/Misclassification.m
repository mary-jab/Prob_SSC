%--------------------------------------------------------------------------
% This function takes the groups resulted from spectral clutsering and the
% ground truth to compute the misclassification rate.
% groups: [grp1,grp2,grp3] for three different forms of Spectral Clustering
% s: ground truth vector
% Missrate: 3x1 vector with misclassification rates of three forms of
% spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------


function [Missrate, FalsePosRate] = Misclassification(groups,s)

n = max(s);

if (n >10)
    Missrate = 1 - compacc(groups, s'); % it is an approximate way to calculate ACC.
else
    Missrate = missclassGroups( groups,s,n ) ./ length(s);
    % it sometimes causes out of memory when nbcluster >10
end

