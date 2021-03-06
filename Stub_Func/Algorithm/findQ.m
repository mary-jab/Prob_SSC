function [QMat] = findQ(clusters, ZKSym, lambda0, lambda1, rho, alpha)
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
QMatThsh = updateQ(ZKSym, clusters, N, lambda0, lambda1, rho, alpha);

QMat = QMatThsh+QMatCls;%(QMatThsh|QMatCls)*(1+eps);%(QMatThsh|QMatCls)*(1+eps);%
end



function [QMatNan] = updateQ(ZKSym, clusters, N, lambda0, lambda1, rho, alpha)
QMatNan = zeros(N,N);
nanIDX = isnan(clusters);
Nan_count = sum(nanIDX);

if Nan_count >0
    NanZ = ZKSym (:,nanIDX);
    
    QMat = findQ_lasso(NanZ,  lambda1, rho, alpha);
    
    
    for i = 1: Nan_count
        QMat(:,i) = QMat(:,i)/max(QMat(:,i));
    end
    QMat(QMat<.3) = 0;
    %     QMat(QMat~=0) = 1;
    
    QMatNan(:,nanIDX) = QMat;
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