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
% iter =1;
% inputOpt.iter = iter;
options.errorPre=clustersErr;
% inputOpt.s = s;
if isfield(options, 'normType')
    normLst = options.normType;
else
    normLst = 2;
end
%% Algorithm
%% first Step
for normType = 2%normLst
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
    for i=1:5
        nameF =strcat('Iter_normType', num2str(normType),'N', num2str(N) );
        if (isfield(options,'sample'))
            nameF =strcat(nameF, 'sample', num2str(options.sample));
        end
        nameF =strcat(nameF, '_iter', num2str(i), '.mat');
        if (exist(fullfile(SAVEPATH,  nameF), 'file'))
            load(fullfile(SAVEPATH,  nameF));
        else
            lambda0_currLst =[.002 .02  .05 .1 .2 .5];%. ]; %.1  .05 0.1 %1;%[10 30 100 500 800 1000];%
            options.errorPre = clustersErr;
            options.itt = i;
            %             inputOpt.GrndTrth = options.GrndTrth;
            missrate(i) = 1;
            lambda0Lst{i}=[];
            lambda1Lst{i}=[];
            thrshPrc{i}=[];
            for lambda0 = lambda0_currLst
                lambda1_currLst = [ lambda0*.0001 ];%, lambda0*.02]; lambda0*.01    30;%
                if i >1
                    lambda1_currLst = [ lambda0*.0001 lambda0*.001 lambda0*.01  lambda0*.01 lambda0*.1];%  lambda0*100];%, lambda0*.02]; lambda0*.01 lambda0*.0001 lambda0*.001
                end %[.001 .01 .1 .6 .9]%[ 1 10 30 80 200 800] ;%
                
                for lambda1 = lambda1_currLst
                    [cZ, cZKSym, cclusters, cclustersErr,CMissrate, cinputOpt] =  mainProcess...
                        (Y, numClass, preQ, preZ, lambda0,lambda1, clusterPre, options, rho, alpha, normType);
                    if CMissrate <= missrate(i)
                        Z=cZ;
                        ZKSym= cZKSym;
                        clusters = cclusters;
                        clustersErr = cclustersErr;
                        inputOpt=cinputOpt;
                        if (CMissrate < missrate(i) && ~isempty(lambda0Lst{i}))
                            lambda0Lst{i} = lambda0;
                            lambda1Lst{i} = lambda1;
                            thrshPrc{i} = inputOpt.thresholdPRC;
                        else
                            lambda0Lst{i}(end+1) = lambda0;
                            lambda1Lst{i}(end+1) = lambda1;
                            thrshPrc{i}(end+1) = inputOpt.thresholdPRC;
                        end
                        missrate(i) = CMissrate;
                    end
                end
            end
            if missrate(i)>=.043
                r =9;
            end
            isNanMat(i) = sum(isnan(clustersErr))/N;
            %               QMat        = findQ_prob(inputOpt.GrndTrth, ZKSym, lambda0, lambda1,  rho, alpha, inputOpt);
            
            QMat        = findQ_prob(clustersErr, ZKSym, lambda0, lambda1,  rho, alpha, inputOpt);
            %             QMat        = findQ(clustersErr, ZKSym, lambda0, lambda1,  rho, alpha);
            %             QMat = findQ_infinityNorm(clustersErr, ZKSym, lambda0, lambda1,  rho, alpha);
            if ( ~isdir(SAVEPATH))
                mkdir(SAVEPATH);
            end
            save(fullfile(SAVEPATH,  nameF), 'Z', 'ZKSym','missrate','isNanMat', 'QMat', 'thrshPrc','lambda0Lst', 'lambda1Lst');
        end
        [sPath] = plotFigure (QMat,ZKSym, missrate(i), options, 11, clustersErr, 0);
        
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
    nameF =strcat('normType', num2str(normType), 'N', num2str(N));
    if (isfield(options,'sample'))
        nameF =strcat(nameF, 'sample', num2str(options.sample));
    end
    nameF =strcat(nameF, '.mat');
    
    
    
    save(fullfile(SAVEPATH,  nameF), 'Z', 'ZKSym','lambda0Lst', 'lambda1Lst', 'thrshPrc', 'missrate', 'QMatLst' );
    mkdir ([ SAVEPATH '/tmp']); movefile([ SAVEPATH '/Iter*.*'], [ SAVEPATH '/tmp/']);
    close all
end
end

function [Z, ZKSym, clusters, clustersErr,missrate, inputOpt] = ...
    mainProcess(Y, numClass, preQ, preZ, lambda0,lambda1, clusterPre, inputOpt, rho, alpha, normType)
[Z, ZKSym] = myLasso(Y, preQ, preZ, lambda0,lambda1, rho, alpha, normType, inputOpt);
[clusters, clustersErr, minErrProb,errorPrbMat,inputOpt] = ...
    Prob_Clustering(Y, ZKSym, numClass, clusterPre, inputOpt);
inputOpt.errorPrbMat=errorPrbMat;
missrate = Misclassification(clusters, inputOpt.GrndTrth); % it is an approximate way to calculate ACC.
end

