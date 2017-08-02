clear all
addpath(genpath('../'))


for k = [4]
    fileName = ['../../../Generate_SyntheticData/genData/SubspaceNoIntersect_affine_Classes_' num2str(k) '_dim_5_ambient_100_error_' num2str(0) '_CVPR_Rene.mat'];
    load (fileName);
    [n, N1, numGrps] = size(Y);
    N= N1*numGrps;
    Y = reshape(( Y(:,:,1:numGrps)),n, N) ;
    
    for errorPrc = [0/100: 10/100:90/100];
        X = Y;
        %% add noise
        sigma =  0.03;
        idx = randperm( length(s), ceil(length(s)*errorPrc));
        for i = 1:length(idx)
            variance = sigma * ((norm(X(:,idx(i)),2)));
            r = (randn(n,1)); r = r-mean(r);%r=r/(std(r));
            gussNoise = ((variance))*r; %Gaussian white noise W
            X(:, idx(i) ) = X(:,idx(i)) +    gussNoise;
        end
        %% run
        opt.SSCrho = 1;
        opt.errorPrc =errorPrc;
        CMatQMat_itt(X, s ,opt);
    end
end

l = 1;