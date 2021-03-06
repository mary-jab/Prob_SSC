if ~isdeployed
    clear all, close all
    addpath(genpath('./'))
addpath(genpath('../back up/S3C-Maryam'))
end

%%  Guassian Noise
DataType{1} = 'syn';
DataType{2} = 'GuassianNoise';
opt.clas = 15;
opt.dim = 5;
opt.amb = 100;
opt.noise = .3;
%% Intersection
% DataType{1} = 'syn';
% DataType{2} = 'intersect';
% opt.clas = 3;
% opt.dim = 10;
% opt.amb = 200;
%% handwritten
% DataType{1} = 'handwritten';

%%
load YaleBCrop025.mat
results_fn =['StrSSC_Faces_tuned_gamma0_results',datestr(now,30),'.mat'];
alpha = 20;
gamma0 =0.3;

opt.ouliers = 1;
opt.affine  = 0;


nSet = [2 3 5 8 10];
for i = 3:length(nSet)
    k = nSet(i);
    idx = Ind{k};
    opt.savePath = ['savedRes/YaleBCrop025/Class_' num2str(k) '/'] ;
    for j = 1:size(idx,1)
        opt.sample = j;
        X = [];
        for p = 1:k
            X = [X Y(:,:,idx(j,p))];
        end
        [D,N] = size(X);
        opt.GrndTrth = s{k}' ;
        opt.k = max(opt.GrndTrth);
        opt.SSCrho = 1;
        %         opt.oulier = 1;
%         load(strcat(opt.savePath, ['normType2N' num2str(length(opt.GrndTrth)) 'sample' num2str(j) '.mat']));
%         if ( (missrate(end)> 0.07))
            Prob_itterative_SSC(X, opt.GrndTrth  ,opt);
%         end
    end
end



% 
% 
% for errorPrc = 00/100: 10/100:60/100
%     opt.noise =errorPrc;
%     for sample = 1:20
%         opt.sample = sample;
%         [Y, opt.GrndTrth , opt.savePath] = ReadData (DataType, opt);
%         k = max(opt.GrndTrth);
%         opt.SSCrho = 1;
%         opt.k = k;
%         
%         
%     end
% end

%         end
%         nameF =strcat('normType', num2str(2), 'N', num2str(size(Y,2)));
%         nameF =strcat(nameF, 'sample', num2str(sample));
%         nameF =strcat(nameF, '.mat');
%         load (fullfile(opt.savePath,  nameF));
%         if (missrate(end)>.05)