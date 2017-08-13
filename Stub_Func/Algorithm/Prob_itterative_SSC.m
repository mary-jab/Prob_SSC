function [result] = Prob_itterative_SSC(Y, s, options)
Debug = 1;

numClass = options.k;
%% run
result = [];
N = length(s);
lambda2= 1.2;
clustersErr = nan(N,1);

errorProb = [];
%% save 1st Step
iter =1;
inputOpt.iter = iter;
inputOpt.errorPre=clustersErr;
inputOpt.s = s;

%% Algorithm
%% first Step
for normType = [2]
    rho = 1.0;
    alpha=1.0;
    lambda0 = .05;%.05;%.05;
    lambda1 = .002*lambda0;
    preQ = ones(N,N);
    clusters = [];
    preZ = zeros(N,N);
    clusterPre = clusters;
    thrshPrc = [];
    
    %% iterative step
    SAVEPATH=strcat(pwd,filesep,options.savePath);
    for i=1:10
        nameF =strcat('Iter_normType', num2str(normType),'N', num2str(N), 'results_iter', num2str(i), '.mat');
        if (exist(fullfile(SAVEPATH,  nameF), 'file'))
            load(fullfile(SAVEPATH,  nameF));
        else
            lambda0_currLst = [.05 .02 0.1]; %.1  .05
            inputOpt.errorPre = clustersErr;
            inputOpt.itt = i;
            inputOpt.GrndTrth = options.GrndTrth;
            missrate(i) = 1;
            lambda0Lst{i}=[];
            lambda1Lst{i}=[];
            thrshPrc{i}=[];
            for lambda0 = lambda0_currLst
                lambda1_currLst = [ lambda0*.0001 ];%, lambda0*.02]; lambda0*.01
                if i >1
                    lambda1_currLst = [ lambda0*.0001 lambda0*.001 lambda0*.01 lambda0*.1];%, lambda0*.02]; lambda0*.01
                end
                
                for lambda1 = lambda1_currLst
                    [cZ, cZKSym, cclusters, cclustersErr,CMissrate, cinputOpt] =  mainProcess...
                        (Y, numClass, preQ, preZ, lambda0,lambda1, clusterPre, inputOpt, rho, alpha, normType);
                    if CMissrate < missrate(i)
                        Z=cZ;
                        ZKSym= cZKSym;
                        clusters = cclusters;
                        clustersErr = cclustersErr;
                        inputOpt=cinputOpt;
                        if (CMissrate < missrate(i) && ~isempty(lambda0Lst{i}))
                            lambda0Lst{i}(end) = lambda0;
                            lambda1Lst{i}(end) = lambda1;
                            thrshPrc{i}(end) = inputOpt.thresholdPRC;
                        else
                            lambda0Lst{i}(end+1) = lambda0;
                            lambda1Lst{i}(end+1) = lambda1;
                            thrshPrc{i}(end+1) = inputOpt.thresholdPRC;
                        end
                        missrate(i) = CMissrate;
                    end
                end
            end
            isNanMat(i) = sum(isnan(clustersErr))/N;
            QMat        = findQ_prob(clustersErr, ZKSym, lambda0, lambda1,  rho, alpha, inputOpt);
%             if sum(isnan(QMat(:)))
%                 yiui=1;
%             end
            %             QMat        = findQ(clustersErr, ZKSym, lambda0, lambda1,  rho, alpha);
%             QMat = findQ_infinityNorm(clustersErr, ZKSym, lambda0, lambda1,  rho, alpha);
            [sPath] = plotFigure (QMat,ZKSym, missrate(i), options, 11, clustersErr, 0);
            if ( ~isdir(SAVEPATH))
                mkdir(SAVEPATH);
            end
            save(fullfile(SAVEPATH,  nameF), 'Z', 'ZKSym','missrate','isNanMat', 'QMat', 'thrshPrc' );
        end
        if (i>2 && isNanMat(i) >= isNanMat(i-1))
            break;
        end
        preQ = QMat;
        preZ = Z;
        clusterPre = clusters;
        QMatLst{i} = QMat;
    end
    
    if ( ~isdir(SAVEPATH))
        mkdir(SAVEPATH);
    end
    nameF =strcat('normType', num2str(normType), 'N', num2str(N),'results.mat');
    save(fullfile(SAVEPATH,  nameF), 'Z', 'ZKSym','lambda0Lst', 'lambda1Lst', 'thrshPrc', 'missrate', 'QMatLst' );
    mkdir ([ SAVEPATH '/tmp']); movefile([ SAVEPATH '/Iter*.*'], [ SAVEPATH '/tmp/']);
    close all
end
end

function [Z, ZKSym, clusters, clustersErr,missrate, inputOpt] = ...
    mainProcess(Y, numClass, preQ, preZ, lambda0,lambda1, clusterPre, inputOpt, rho, alpha, normType)
[Z, ZKSym] = myLasso(Y, preQ, preZ, lambda0,lambda1, rho, alpha, normType);
[clusters, clustersErr, minErrProb,errorPrbMat,inputOpt] = ...
    Prob_Clustering(Y, ZKSym, numClass, clusterPre, inputOpt);
inputOpt.errorPrbMat=errorPrbMat;
missrate = Misclassification(clusters, inputOpt.GrndTrth); % it is an approximate way to calculate ACC.
end

