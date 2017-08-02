function [missrate,normLst, lambdaBstLst, grps] = CMatQMat_itt(Y, clusterIDX ,opt)
%% Initialization
[n,N] = size(Y);
rho =opt.SSCrho; %rho = 1 be default;
nbcluster = max(clusterIDX);
QMat = ones(N, N);
QMat_old=1-QMat;
lambda00 = 1;
patternQ = [];
for i=1:N,
    for j=1:N,
        if clusterIDX(i)==clusterIDX(j)
            patternQ(i,j) =1+eps;
        end
    end
end

%% %%%%%%%%%%%%%%%%%
missrate=[]; normLst=[];
lambda2 = .5;
lambda0Max = .7;%*(1/computeLambda_mat(Y));
lambda0 = .05;
for iter=1:10
    t_start = tic;
    
    if iter ==1;
        lstLambda0 = .05; lstLambda1=0;
        lstLambda2 =0;
    else
        lstLambda0 =[.05 : .1 : lambda0Max];
        lambda0 = lstLambda0(1);
        lstLambda1 = [lambda0/2:lambda0/2:lambda0*5];
        lstLambda2 =1.1:.2:2.5;
    end
    %% solve the re-weighted SSC problem
    missrateT = []; normLstT= []; lambda0T = []; lambda1T = []; lambda2QT =[];ZKSym_tmp = [];QMat_tmp = [];
    for ii = 1:length(lstLambda0)
        lambda0 = lstLambda0(ii);
        for jj = 1:length(lstLambda2)
            if (iter>1)
                lambda2 = lstLambda2(jj);
                for i = 1 : N,
                    [~, id] = sort(ZKSym_old (:,i), 'descend' );
                    QMat(:,i) = zeros(N,1);
                    QMat(id(1:2),i) = 1+eps;
                    QMat(:,i) = (QMat(:,i)|(ZKSym_old (:,i) >  lambda2 *std(ZKSym_old(ZKSym_old(:,i)>0,i))))*(1+eps);
                end
            end
            for kk = 1:length(lstLambda1)
                lambda1 = lstLambda1(kk);
                %% calculate classification error rate
                [ZKSym,CAbs, missrateT(end+1)] =  findCoefCluster (Y, QMat, clusterIDX, lambda0, lambda1);
                normLstT(end+1) = norm(~patternQ.*ZKSym,1);
                
                lambda0T(end+1)  = lambda0;
                lambda1T(end+1)  = lambda1;
                lambda2QT(end+1)  = lambda2;
                ZKSym_tmp{end+1} = ZKSym;
                QMat_tmp{end+1}  = QMat;
            end
        end
    end
    [v,id]           = min(missrateT); ids = find(missrateT==v);
    missrate(iter)   = missrateT(id);
    normLst{iter}    = normLstT(ids);
    lambdaBst0{iter} = lambda0T(ids);
    lambdaBst1{iter} = lambda1T(ids);
    lambdaBst2Q{iter} = lambda2QT(ids);
    QMat = QMat_tmp{id};
    
    %% threshold on Z
    toc(t_start);
    if max(max(abs(QMat-QMat_old)))<eps
        save(['saveData/Details/Results_Prc_' num2str(opt.errorPrc) '_class_' num2str(nbcluster) '.mat'],...
            'iter', 'missrate', 'normLst', 'lambdaBst0', 'lambdaBst1', 'lambdaBst2Q');
        break;
    end
    QMat_old = QMat_tmp{id};
    ZKSym_old = ZKSym_tmp{(id)};
    
end

save(['saveData/Details/Results_Prc_' num2str(opt.errorPrc) '_class_' num2str(nbcluster) '.mat'],...
    'iter', 'missrate', 'normLst', 'lambdaBst0', 'lambdaBst1', 'lambdaBst2Q');
end

