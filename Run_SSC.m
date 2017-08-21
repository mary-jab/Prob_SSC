if ~isdeployed
    clear all, close all
    addpath(genpath('./'))
end

%%  Guassian Noise
DataType{1} = 'syn';
DataType{2} = 'GuassianNoise';
opt.clas = 15;
opt.dim = 5;
opt.amb = 100;
opt.noise = .3;
%% Intersection
DataType{2} = 'intersect';
opt.clas = 2;
opt.dim = 10;
opt.amb = 200;
opt.ouliers = 0;
opt.affine = 0;


%% handwritten
% DataType{1} = 'handwritten';

%%
normType = 2;
N = opt.clas *100;
cnt = 0;
for cls = 2:2
    opt.clas = cls;
    for errorPrc = 40/100: 10/100:60/100
        cnt = cnt+1;
        opt.noise =errorPrc;
        
        loadPath = strcat('savedRes\Intersect\cls' ,num2str(cls) ,'\Subspace_noise_', num2str(errorPrc), '\ambiant200');
        for sample = 1:50
            nameF =strcat('normType', num2str(normType), 'N', num2str(N));
            nameF =strcat(nameF, 'sample', num2str(sample));
            nameF =strcat(nameF, '.mat');
            load (fullfile(loadPath,  nameF));
            misArr(sample, :) = [missrate, ones(1, 10-length(missrate))* missrate(end)];
            misEnd(sample) = missrate(end);
        end
        [a,  idx] = sort(misEnd);
        misEnd(idx(21:end)) = [];
        misArr((idx(21:end)),:) = [];
        idxList(cnt, :) = idx(1:20);
        
       
        for sample =idxList(cnt,18:20)
            opt.sample  = idxList(1);
            [Y, opt.GrndTrth , opt.savePath] = ReadData (DataType, opt);
            k = max(opt.GrndTrth);
            opt.SSCrho = 1;
            opt.k = k;
                        opt.sample = sample;

                        opt.savePath = strcat(opt.savePath, '/SSC');
            itterative_SSC(Y, opt.GrndTrth  ,opt);
        end
    end
end

