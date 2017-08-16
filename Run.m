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
% DataType{1} = 'syn';
% DataType{2} = 'intersect';
% opt.clas = 3;
% opt.dim = 10;
% opt.amb = 200;
%% handwritten
% DataType{1} = 'handwritten';

%%

for errorPrc = 00/100: 10/100:60/100
    opt.noise =errorPrc;
    for sample = 1:20
        opt.sample = sample;
        [Y, opt.GrndTrth , opt.savePath] = ReadData (DataType, opt);
        k = max(opt.GrndTrth);
        opt.SSCrho = 1;
        opt.k = k;
        
        Prob_itterative_SSC(Y, opt.GrndTrth  ,opt);
        
    end
end
        
                %         end
        %         nameF =strcat('normType', num2str(2), 'N', num2str(size(Y,2)));
        %         nameF =strcat(nameF, 'sample', num2str(sample));
        %         nameF =strcat(nameF, '.mat');
        %         load (fullfile(opt.savePath,  nameF));
        %         if (missrate(end)>.05)