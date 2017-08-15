
noise = .8;

loadPath = strcat(pwd,'\savedRes\GuassianNois_Cls_15\Subspace_noise_', num2str(noise), '\ambiant100');
normType = 2;
N = 150;
misArr = zeros(20,10);
for sample = 1:10
    nameF =strcat('normType', num2str(normType), 'N', num2str(N));
        nameF =strcat(nameF, 'sample', num2str(sample));
    nameF =strcat(nameF, '.mat');
    load (fullfile(loadPath,  nameF));
    misArr(sample, :) = [missrate, ones(1, 10-length(missrate))* missrate(end)];
    misEnd(sample) = missrate(end);
end

mean(misEnd)

figure, hold on 
c = 'rgbmkcrgbmkcrgbmkc';
for sample = 1:5
    plot (misArr(sample, :), 'color', c(sample) )
end