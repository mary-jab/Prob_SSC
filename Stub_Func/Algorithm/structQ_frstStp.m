
function [QMat] = structQ_frstStp(clusters)
N = length(clusters);
QMat = zeros(N,N);
%% clusters
for i=1:N,
    for j=1:N,
        if clusters(i)==clusters(j)
            QMat(i,j) =1+eps;
        end
    end
end
end