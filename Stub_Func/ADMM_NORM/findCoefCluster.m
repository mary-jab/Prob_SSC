function [CKSym,CAbs, missrate] =  findCoefCluster (Y, Q, classIdx, lambda0, lambda1)

C = lasso_Q(Y, Q, lambda0, lambda1, 1.0, 1.0);
[CKSym,CAbs] = BuildAdjacency(C);
%% calculate classification error rate
grps = SpectralClustering(CKSym, max(classIdx));
if xor(iscolumn(grps), iscolumn(classIdx))
    classIdx = classIdx';
end
%                 missrateT(end+1) = 1 - compacc(grps, idx); % it is an approximate way to calculate ACC.
missrate = Misclassification(grps, classIdx); % it sometimes causes out of memory when nbcluster >10

end



% C = lasso_Q(Y, Q, .4, lambda1, 1.0, 1.0);
% [CKSym,CAbs] = BuildAdjacency(C);