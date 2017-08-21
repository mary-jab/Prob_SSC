clear all, close all
saveArr = []; cnt = 0;
for noise = 40/100: 10/100:60/100
    % loadPath = strcat(pwd,'\savedRes\GuassianNois_Cls_15\Subspace_noise_', num2str(noise), '\ambiant100');
    % N = 150;
    cnt = cnt+1;
    
    cls = 4;
    loadPath = strcat(pwd,'\savedRes\Intersect\cls' ,num2str(cls) ,'\Subspace_noise_', num2str(noise), '\ambiant200\SSC\');
    N = cls*100;
    % cls = 2;
    % loadPath = strcat(pwd,'\savedRes\\YaleBCrop025\Class_', num2str(cls));
    % N = cls * 64;
    
    normType = 2;
    
    dataLstName=[];
    imgDir = dir (loadPath);
    for i=1:length(imgDir)
        if ~(imgDir(i).isdir)
            dataLstName{end+1} = [imgDir(i).name];
        end
    end
    
    
    % misArr = zeros(20,10);
    for sample = 1:length(dataLstName)
        nameF =dataLstName{sample};
          load (fullfile(loadPath,  nameF));
        misArr(sample, :) = [missrate, ones(1, 10-length(missrate))* missrate(end)];
        misEnd(sample) = missrate(end);
    end
 
    
    
    saveArr(cnt, 1) = mean(misEnd);
    SSC(cnt,1) = mean(misArr(:,1));
    % saveArr(cnt, 2) = median(misEnd);
    
    % figure, hold on
    % c = 'rgbmkcrgbmkcrgbmkcrgbmkcrgbmkcrgbmkc';
    % lst = [];
    % for sample = 1:20
    %     if (misEnd(sample)>0.3)
    %         lst(end+1) = sample;
    %     plot (misArr(sample, :), 'color', c(sample) )
    %     end
    % end
end
k = 1;