function clusters= NCutCluster(A, maxCluster)
N = size(A,1);
clusters = zeros(N, 1);

[NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(A, maxCluster);
for i=1:N
    [~,clusters(i)] = find(NcutDiscrete(i,:)==1);
end