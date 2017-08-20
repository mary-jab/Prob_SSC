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
list= [ 44    50    17    33    12    18    43     3    37     9    16     4    31    29    30    45    23    24     5    35;
        35    49    18    33     9    25    46     4    37    31    44    14    12     2    11    34    45    27    22    38];
cnt = 0;
for cls = 4:4
    opt.clas = cls;
    for errorPrc = 50/100: 10/100:60/100
            cnt = cnt+1;
        opt.noise =errorPrc;
        for sample =list(cnt,:)
            opt.sample = sample;
            [Y, opt.GrndTrth , opt.savePath] = ReadData (DataType, opt);
            k = max(opt.GrndTrth);
            opt.SSCrho = 1;
            opt.k = k;
            opt.savePath = strcat(opt.savePath, '/SSC');
            itterative_SSC(Y, opt.GrndTrth  ,opt);
        end
    end
end

%         end
%         nameF =strcat('normType', num2str(2), 'N', num2str(size(Y,2)));
%         nameF =strcat(nameF, 'sample', num2str(sample));
%         nameF =strcat(nameF, '.mat');
%         load (fullfile(opt.savePath,  nameF));
%         if (missrate(end)>.05)