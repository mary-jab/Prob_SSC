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
for normType = [2, 1]
    rho = 1.0;
    alpha=1.0;
    lambda0 = .05;%.05;%.05;
    lambda1 = .002*lambda0;
    preQ = ones(N,N);
    clusters = [];
    preZ = zeros(N,N);
    clusterPre = clusters;
    
    
    %% iterative step
    for i=1:10
        lambda0_currLst = .1;
        
        inputOpt.errorPre = clustersErr;
        inputOpt.itt = i;
        inputOpt.GrndTrth = options.GrndTrth;
        
        missrate(i) = 1;
        lambda0Lst{i}=[];
        lambda1Lst{i}=[];
        thrshPrc{i}=[];
        for lambda0 = lambda0_currLst
            lambda1_currLst = [ lambda0*.0001, lambda0*.0005 lambda0*.001 lambda0*.01];%, lambda0*.02];
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
        QMat        = findQ(clustersErr, ZKSym, lambda0, lambda1,  rho, alpha);
        %         [sPath] = plotFigure (QMat,ZKSym, missrate(i), options, 11, clustersErr, 0);
        
        if (i>1 && isNanMat(i) >= isNanMat(i-1))
            break;
        end
        preQ = QMat;
        preZ = Z;
        clusterPre = clusters;
        QMatLst{i} = QMat;
        
    end
    if ~isdeployed
        mkdir(options.savePath)
        save ([options.savePath 'normType' num2str(normType) 'results.mat'],...
            'Z', 'ZKSym','lambda0Lst', 'lambda1Lst', 'thrshPrc', 'missrate', 'QMatLst' );
    else
        SAVEPATH=strcat(pwd,filesep,'output');
        if ( ~isdir(SAVEPATH))
            mkdir(SAVEPATH);
        end
        nameF =strcat('normType', num2str(normType), 'results.mat');
        save(fullfile(SAVEPATH,  nameF), 'Z', 'ZKSym','lambda0Lst', 'lambda1Lst', 'thrshPrc', 'missrate', 'QMatLst' );
        %         nameF = ['normType' num2str(normType) 'results.mat'];
        %         save(fullfile(ctfroot, mfilename, nameF), 'Z', 'ZKSym','lambda0Lst', 'lambda1Lst', 'thrshPrc', 'missrate', 'QMatLst' , '-append');
    end
    close all
end
end

function [Z, ZKSym, clusters, clustersErr,missrate, inputOpt] = ...
    mainProcess(Y, numClass, preQ, preZ, lambda0,lambda1, clusterPre, inputOpt, rho, alpha, normType)
[Z, ZKSym] = myLasso(Y, preQ, preZ, lambda0,lambda1, rho, alpha, normType);
[clusters, clustersErr, minErrProb, errorPrbMat,inputOpt] = ...
    Prob_Clustering(Y, ZKSym, numClass, clusterPre, inputOpt);
missrate = Misclassification(clusters, inputOpt.GrndTrth); % it is an approximate way to calculate ACC.
end



%
%
%
%
% [missrateT, missrateConf] = MissrateMat(clustersList,  options, errorProb{end});
%
%
%
%
%
% inputOpt.initial = 1;
% if numClass>2
%     [clustersList,  clusErrList,  errScr, errorProb{end+1}, errTrsh]...
%         = clusterDS_MaryamNcut_3Cls(options.Y, options.CKSym, numClass, preClusters, inputOpt);
% else
%     [clustersList,  clusErrList,  errScr, errorProb{end+1}, errTrsh]...
%         = clusterDS_MaryamNcut(options.Y, options.CKSym, numClass, preClusters, inputOpt);
% end
% % [clustersList,  clusErrList,  errScr, errorProb{end+1}, errTrsh] = clusterDS_MaryamNcut(options.Y, options.CKSym, numClass, preClusters, inputOpt);
% % missrateT   = Misclassification(clustersList(:,end), s); % it is an approximate way to calculate ACC.
% [missrateT, missrateConf] = MissrateMat(clustersList,  options, errorProb{end});
% % missrateT = missrateConf(1);
%
% %%%%%%%%%%%%%%%%%%%% delete !!!!!!!!!!!!!!!!!
% % clusErrList = clustersList;
% % clusErrList (options.GrndTrth~=options.confusionGrndTrth | options.confusionErrProb>.45)= NaN;
% %%
% QMat        = findQ(clusErrList, options.CKSym, lambda2);
%
% %%
% errorPre = clusErrList;
% result.misClassification(iter) = options.missrate;%missrateT;
% result.misCls2(iter)           = missrateConf(1);
% result.misClsNaN(iter)         = missrateConf(2);
% result.AccPR(iter)             = missrateConf(3);
%
% result.ClusterLst(:,iter)      = clustersList;
% result.clsErrLst(:,iter)       = clusErrList;
% result.errScrLst(:,iter)       = errScr;
% result.errProb{iter}           = errorProb;
% result.ERRlst(iter)            = estERR(options.CMat, Y);
% result.featureErr(iter)        = estERR_feature(options.CKSym, clustersList);
% result.QMat{iter}              = QMat;
% result.ZKSym{iter}             = options.CKSym;
% result.errTrsh(iter)           = errTrsh;
% bestClas = result;
% % if Debug
% [sPath] = plotFigure (QMat, result.misClassification(iter), options, iter, errorPre);
% % end
% save ([sPath '/detail_' num2str(options.k) '_Err_' num2str(options.errorPrc) '_smpl_' ...
%     num2str(options.itt) '_itr_' num2str(iter) '.mat' ], 'bestClas');
%
% Counter = 1;
% bestClas=[];
% inputOpt.initial=0;
% %% 1st step
% %% Initialization
% preZ = options.CMat; %zeros(size(D,2));
%
% swing = 0;
%
% while (iter<options.numOfIteration)
%
%     if iter ==1
%         lambda = [ .03 .06 .1 .2 .3 .4];%.0002;
%     else
%         lambda =[.06 .1 .2];%[ .02 .06 .1 .2 .3 .4] ;% [ .0002  .02]%[.0002 .0008];
%     end
%
%     %% next step
%     iter = iter+1;
%     for t = 1:2
%         for lambda0 = lambda; %.02 .2
%             if iter>1   %% find Q
%                 lambda1List = [lambda0, 10*lambda0 ];%[2*lambda0, 5*lambda0:5*lambda0:30*lambda0];%lambda0/10:lambda0:30*lambda0;%10*lambda0;%
%             else
%                 lambda1List = 0;
%             end
%             for lambda1 = lambda1List
%                 errorProb = [];
%                 %% find C
%                 rho = 1.0;
%                 alpha=1.0;
%                 [Z, ZKSym] = myLasse(Y, QMat,preZ,lambda0,lambda1, rho, alpha);
%                 %% dominant Set
%
%                 inputOpt.errorPre=errorPre;
%                 inputOpt.iter = iter;
%                 if numClass>2
%                     [clustersList,  clusErrList,  errScr, errorProb{end+1}, errTrsh]...
%                         = clusterDS_MaryamNcut_3Cls(Y, ZKSym, numClass, preClusters, inputOpt);
%                 else
%                     [clustersList,  clusErrList,  errScr, errorProb{end+1}, errTrsh]...
%                         = clusterDS_MaryamNcut(Y, ZKSym, numClass, preClusters, inputOpt);
%                 end
%                 %% miss class
%                 %         missrateT = 1 - compacc(clustersList(:,end), s); % it is an approximate way to calculate ACC.
%                 %         missrateT = Misclassification(clustersList(:,end), s); % it sometimes causes out of memory when nbcluster >10
%                 %missrateT = Misclassification(clustersList(:,end),s); % it is an approximate way to calculate ACC.
%
%                 [missrateT, missrateConf] = MissrateMat(clustersList(:,end), options,errorProb{end});
%                 %         missrateT = missrateConf(1);
%
%
%                 IDX = 1:length(missrateT);
%                 minID       = IDX(missrateT==min(missrateT));
%                 clusters    = (clustersList(:,minID(1)));
%                 clsErr      = clusErrList(:, minID(1));
%                 bestClas.lambda0(Counter)       = lambda0;
%                 bestClas.lambda1(Counter)       = lambda1;
%                 bestClas.misClas(Counter)       = missrateT(minID(1));
%
%                 bestClas.misClas2(Counter)      = missrateConf(1);
%
%                 bestClas.misClasNaN(Counter)    = missrateConf(2);
%                 bestClas.AccPR(Counter)         = missrateConf(3);
%
%                 bestClas.clusters(:,Counter)    = clusters(:,1);
%                 bestClas.clusErrList(:,Counter) = clsErr(:,1);
%                 bestClas.errScr(:,Counter)      = errScr(:, minID(1));
%                 bestClas.ERR(Counter)           = estERR(Z, Y);
%                 bestClas.featureErr(Counter)    = estERR_feature(ZKSym, clusters);
%                 bestClas.errProb{Counter}       = errorProb{ minID(1)};
%                 bestClas.errTrsh(Counter)       = errTrsh;
%                 ZKList{Counter}                 = ZKSym;
%                 ZList{Counter}                  = Z;
%                 Counter = Counter+1;
%
%             end
%         end
%         IDX = 1:Counter-1;
%         minMisCls   = IDX(min(bestClas.misClas )==bestClas.misClas);
%         minID       = minMisCls(min(bestClas.misClasNaN(minMisCls))==bestClas.misClasNaN(minMisCls));
%         ZKSymPre    = ZKList{minID(end)};
%         preZ        = ZList{minID(end)};
%         newClusters = bestClas.clusters(:,minID(end));
%         errorPre    = bestClas.clusErrList(:,minID(end));
%     end
%
%     result.misClassification(iter) = bestClas.misClas(minID(end));
%     result.misCls2(iter)           = bestClas.misClas2(minID(end));
%     result.misClsNaN(iter)         = bestClas.misClasNaN(minID(end));
%     result.AccPR(iter)             = bestClas.AccPR(minID(end));
%     result.ClusterLst(:,iter)      = newClusters;
%     result.clsErrLst(:,iter)       = bestClas.clusErrList(:,minID(end));
%     result.errScrLst(:,iter)       = bestClas.errScr(:,minID(end));
%     result.errProb{iter}           = bestClas.errProb{minID};
%     result.lambda0List(iter)       = bestClas.lambda0(minID(end));
%     result.lambda1Iter(iter)       = bestClas.lambda1(minID(end));
%     result.ERRlst(iter)            = bestClas.ERR(minID(end));
%     result.featureErr(iter)        = bestClas.featureErr(minID(end));
%     result.QMat{iter}              = QMat;
%     result.ZKSym{iter}             = ZKSymPre;
%     result.errTrsh(iter)           = bestClas.errTrsh(minID(end));
%     QMatPre      = QMat;
%
%     %     errorPre = options.GrndTrth;
%     %     errorPre (options.GrndTrth~=options.confusionGrndTrth)= NaN;
%     %         errorPre (newClusters~=options.confusionGrndTrth)= NaN;
%
%
%     [QMat]       = findQ(errorPre, ZKSymPre, lambda2);
%     %% save
%     if Debug
%         [sPath] = plotFigure (QMat, result.misClassification(iter), options, iter, errorPre);
%     end
%     save ([sPath '/detail_' num2str(options.k) '_Err_' num2str(options.errorPrc) '_smpl_' ...
%         num2str(options.itt) '_itr_' num2str(iter) '.mat' ], 'bestClas');
%
%     %% next step
%     bestClas = []; Counter = 1;
%     if ( sum(isnan(result.clsErrLst(:,iter))) >= sum(isnan(result.clsErrLst(:,iter-1))))
%         swing = swing+2;
%     end
%
%     if      ...%(Misclassification(preClusters, prePreClusters)<eps ||...
%             (Misclassification(preClusters, newClusters)<eps)...
%             || swing>1%( sum(isnan(result.clsErrLst(:,iter))) > sum(isnan(result.clsErrLst(:,iter-1)))) ...
%         %              ||  Misclassification(result.clsErrLst(:,iter) , result.clsErrLst(:,iter-1)) <eps
%         %  ||  max(max(abs(QMat-QMatPre)))<1
%         break;
%     end
%     prePreClusters = preClusters;
%     preClusters = newClusters;
% end
%
% if ~any(clusters)
%     clusters= fillClusters(Y, clusters);
% end
%
%
%
% end
%
%
% %% stub Func
%
%
%
%
%
% %% calculate estimation error
% function ERR = estERR(Z, Y)
% D = Y * Z;
% ERR = norm(D-Y, 'fro')/ norm(Y, 'fro');
% end
%
%
% %%
% function [CSym, normC]=NormalizeC(CMat)
%
% N = size(CMat,1);
% CMat = abs(CMat);
% CMat_Sync = CMat+CMat';
% for i=1:N
%     normC(:,i) = CMat(:,i)/(max(CMat(:,i)));
%     CSym(:,i) = CMat_Sync(:,i)/(max(CMat_Sync(:,i))+eps);
% end
% end
%
