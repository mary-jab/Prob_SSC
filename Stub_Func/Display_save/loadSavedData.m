clear all
for noise =  50/100: 10/100:60/100

% loadPath = strcat(pwd,'\savedRes\GuassianNois_Cls_15\Subspace_noise_', num2str(noise), '\ambiant100');
% N = 150;


cls = 2;
loadPath = strcat(pwd,'\savedRes\Intersect\cls' ,num2str(cls) ,'\Subspace_noise_', num2str(noise), '\ambiant200');
N = cls*100;
% cls = 2;
% loadPath = strcat(pwd,'\savedRes\\YaleBCrop025\Class_', num2str(cls));
% N = cls * 64;

normType = 2;

% misArr = zeros(20,10);
for sample = 1:20
    nameF =strcat('normType', num2str(normType), 'N', num2str(N));
        nameF =strcat(nameF, 'sample', num2str(sample));
    nameF =strcat(nameF, '.mat');
    load (fullfile(loadPath,  nameF));
    misArr(sample, :) = [missrate, ones(1, 10-length(missrate))* missrate(end)];
    misEnd(sample) = missrate(end);
end


mean(misEnd)

figure, hold on 
c = 'rgbmkcrgbmkcrgbmkcrgbmkcrgbmkcrgbmkc';
lst = [];
for sample = 1:20
    if (misEnd(sample)>0.3)
        lst(end+1) = sample;
    plot (misArr(sample, :), 'color', c(sample) )
    end
end
end