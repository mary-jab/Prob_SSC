function [clusters, clustersErr, minErrProb, errorPrbMat, options]...
    = Prob_Clustering(Y, ZKSym, maxCluster, preClusters, options)
% computThresh = 0;
% if nargin < 6
%     computThresh = 1;
% %     theshold = 1/(maxCluster*maxCluster)/(options.itt);%*options.itt); % maxCluster;.8;%
% end
N = size(ZKSym,1);
%clustering
clusters= NCutCluster(ZKSym, maxCluster);
%clusters = SpectralClustering(A, maxCluster, preClusters);

%% projection error
% [errorPrbMat, errRatio, errRatioIDX, minErrProb, pickedClstr]= calProjErr (Y, clusters, maxCluster);
[errorPrbMat, errRatio, errRatioIDX, minErrProb, pickedClstr]= CalcNormErr(ZKSym, clusters, maxCluster);

i = 1;
minMis = 1;
for prc = [.95:-.01:.05]
    % if isfield(options, 'BaseTheshold')==0
    sortedArr = sort(errRatio);
    options.BaseTheshold  = sortedArr(uint16(prc*N));
    % end
    theshold =options.BaseTheshold/(options.itt);%sortedArr(uint16((.9-(options.itt-1)*.1)*N));% *maxCluster);
    
    diffIdx1 = (pickedClstr~=clusters);
    diffIdx2 = errRatio<theshold ;%;| minErrProb >= 1/(1.5*maxCluster);
%     sum(diffIdx2)
    clustersErr= clusters ;
    clustersErr( isnan(options.errorPre) & diffIdx2) = NaN;
    % if (options.itt<5)
    clustersErr(diffIdx1) = NaN;
    
    id = ~isnan(clustersErr);
    missrate(i) = Misclassification(clustersErr(id), options.GrndTrth(id)); % it is an approximate way to calculate ACC.
    if minMis>=missrate(i);
        minPrc = prc;
        minMis = missrate(i);
        minCls = clustersErr;
    end
    i = i+1;
end
clustersErr = minCls;
options.thresholdPRC = minPrc;

% clusters = pickedClstr;
% end

end
function [errorPrbMat, errRatio, errRatioIDX, minErrProb, pickedClstr]= CalcNormErr(ZKSym, clusters, maxCluster)
errRatioIDX = [];
N = length(clusters);

errorPrbMat = zeros(N, maxCluster);
pickedClstr = zeros(N,1);
minErrProb =  zeros(N,1);
errRatio =  zeros(N,1);
for i=1:N
    currentCls = clusters;
    currentCls(i) = 0;
    D = calNormRatio(ZKSym(:,i), currentCls, maxCluster);
    [~,C] = max(D);
    id = true(1,maxCluster); id(C) = 0;
    errRatio(i) = min(abs(D(id)-D(~id)));
    minErrProb(i) = 1-D(~id);
    pickedClstr(i)= C;
    errorPrbMat(i,:) = D;
end
end

function [D] = calNormRatio(ZKSym_row, currentCls, maxCls)
normZi = norm(ZKSym_row,1);
for j=1:maxCls
    idx= j==currentCls;
    CurrentA = ZKSym_row(idx);
    D(j) = norm(CurrentA,1)/normZi;
    %errorPrbMat(i,j) = 1- norm(CurrentA,1)/norm(A(:,i),1);
end
end

