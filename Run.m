clear all, close all
addpath(genpath('./'))


%%  Guassian Noise
DataType{1} = 'syn';
DataType{2} = 'GuassianNoise';
opt.clas = 15;
opt.dim = 5;
opt.amb = 100;
% opt.noise = .2;
%% Intersection
% DataType{1} = 'syn';
% DataType{2} = 'intersect';
% opt.clas = 3;
% opt.dim = 10;
% opt.amb = 200;
%% handwritten
DataType{1} = 'handwritten';

%%
opt.samNum = 1;

for errorPrc = 20/100: 10/100:90/100;
    opt.noise = errorPrc;
    [Y, runOpt.GrndTrth , runOpt.savePath] = ReadData (DataType, opt);
    opt.clas = max(runOpt.GrndTrth);
    runOpt.SSCrho = 1;
    runOpt.errorPrc =errorPrc;
    runOpt.k = opt.clas;
    
    Prob_itterative_SSC(Y, runOpt.GrndTrth  ,runOpt);
end
