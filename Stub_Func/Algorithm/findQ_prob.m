function [QMat] = findQ_prob(clusters, ZKSym, lambda0, lambda1, rho, alpha, inputOpt)
N = length(clusters);
QMatCls = zeros(N,N);
%% clusters
for i=1:N,
    for j=1:N,
        if clusters(i)==clusters(j)
            QMatCls(i,j) =1+eps;
            QMatCls(j,i) =1+eps;
        end
    end
end
%% threshold
QMatThsh = updateQ(ZKSym, clusters, N, inputOpt);

QMat = QMatThsh+QMatCls;%(QMatThsh|QMatCls)*(1+eps);%(QMatThsh|QMatCls)*(1+eps);%
end

function [QMatNan] = updateQ2(ZKSym, clusters, N, inputOpt)
maxCluster = max(clusters);
QMatNan = zeros(N,N);
nanIDX = find(isnan(clusters));
Nan_count = sum(nanIDX);
errorPrbMat = inputOpt.errorPrbMat;
for m=1:length(nanIDX)
    i = nanIDX(m);
    
    for j = 1:maxCluster
        curIdx = clusters==j;
        QMatNan(i,curIdx) = inputOpt.errorPrbMat(i,j);
        QMatNan(curIdx, i) = inputOpt.errorPrbMat(i,j);
    end
end

end

function [QMatNan] = updateQ(ZKSym, clusters, N, inputOpt)
maxCluster = max(clusters);
QMatNan = zeros(N,N);
nanIDX = find(isnan(clusters));
Nan_count = sum(nanIDX);
errorPrbMat = inputOpt.errorPrbMat;
for m=1:length(nanIDX)
    i = nanIDX(m);
    [~, firstMaxId] = max(errorPrbMat(i,:));
    mat = firstMaxId;
    for it=1:2
        if maxCluster>2
            errorPrbMat(i,mat(end)) = -1;
            if  ~isnan(max(errorPrbMat(i,:)))
                [~, ScndMaxId] = max(errorPrbMat(i,:));
                mat = [mat, ScndMaxId] ;
            end
        end
    end
    if (sum(inputOpt.errorPrbMat(i,mat))~=0)
        a(mat) = inputOpt.errorPrbMat(i,mat)/sum(inputOpt.errorPrbMat(i,mat));
        for j = mat
            curIdx = clusters==j;
            QMatNan(i,curIdx) = a(j);%inputOpt.errorPrbMat(i,j);
            QMatNan(curIdx, i) = a(j);%inputOpt.errorPrbMat(i,j);
        end
    end
end

end



function QMat = updateCVX_itt(NanZ, lambda0,lambda1 ,N, Nan_count)

for i = 1:Nan_count
    cvx_begin quiet;
    variable q(N)
    expressions YY
    YY=  norm( (1-q) .* NanZ(:,i), 2);
    minimize( sum( (YY))+ lambda0* norm(q,1));%+ lambda2 *sum(X) );%) sum(X(:))
    subject to
    %sum( QMat,1)==1;
    q>=0;
    q<=1;
    cvx_end;
    QMat(:,i) = q;
    
end

end

function QMat = updateCVX(NanZ, lambda0,lambda1 ,N, Nan_count)


cvx_begin quiet;
variable QMat(N,Nan_count)
expressions YY(Nan_count)
for i = 1 : Nan_count,
    YY(i)= lambda1 * norm( (1-QMat(:,i)) .* NanZ(:,i), 2);
end
minimize( sum( (YY))+ lambda0* norm(QMat,1));%+ lambda2 *sum(X) );%) sum(X(:))
subject to
%sum( QMat,1)==1;
QMat>=0;
QMat<=1;
cvx_end;

end