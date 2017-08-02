function [C, CKSym, clusters, missrate, QMat] = firstStepSSC(Y,GrndTrth, lambda0, rho, alpha)
numClass = max(GrndTrth);
[C, ~] = lasso_itt(Y, lambda0, rho, alpha);

[CKSym,CAbs] =  BuildAdjacency(thrC(C,1));
clusters = SpectralClustering(CKSym, numClass);
missrate = Misclassification(clusters, GrndTrth); % it is an approximate way to calculate ACC.
QMat = structQ_frstStp(clusters);

end



