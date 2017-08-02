function [errorPrbMat, errRatio, errRatioIDX, minErrProb, pickedClstr]=...
    calProjErr (Y, clusters, maxCluster)

N = length(clusters);

errorPrbMat = zeros(N, maxCluster);
pickedClstr = zeros(N,1);
minErrProb =  zeros(N,1);
errRatio =  zeros(N,1);
parfor i=1:N
    idxLst = true(N,1);
    idxLst(i) = 0;
    kappa = 1e-6;% parameter for CRC, which could be used for avioding over-fitting
    Proj_M = (Y(:,idxLst)'*Y(:,idxLst)+kappa*eye(size(Y(:,idxLst),2)))\Y(:,idxLst)';
    [pickedClstr(i),~, E] = CRC_IDcheck(Y(:,idxLst),Proj_M,Y(:,~idxLst),clusters(idxLst), maxCluster); 
    sumE = sum(E(~isnan(E) & ~isinf(E)));
    errProb = E/sumE;
    errorPrbMat(i,:) = errProb;
    
    [~, idx] = sort(E);
    errRatio(i) = errProb(idx(2))-errProb(idx(1));
    errRatioIDX(i,:) = [idx(1) idx(2)];
    minErrProb(i) = E(clusters(i))/sumE;
    
end

end