clear all, close all
saveArr = []; cnt = 0; S3C_list = [];
for cls = 2:2
    for noise = 40/100: 10/100:90/100
        % loadPath = strcat(pwd,'\savedRes\GuassianNois\Cls_15\Subspace_noise_', num2str(noise), '\ambiant100');
        % N = 150;
        cnt = cnt+1;
        loadPath = strcat(pwd,'\savedRes\Intersect\cls' ,num2str(cls) ,'\Subspace_noise_', num2str(noise), '\ambiant200');
        N = cls*100;
        % cls = 2;
        % loadPath = strcat(pwd,'\savedRes\\YaleBCrop025\Class_', num2str(cls));
        % N = cls * 64;
        normType = 2;
        
        % misArr = zeros(20,10);
        for sample = 1:30
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
        
        probSCC_list(cnt, 1) = mean(misEnd);
        SCC_list(cnt,1) = mean(misArr(:,1));
        
        misS3C = [];
        for sample = idxList(cnt,:)
            nameF =strcat('normType', num2str(normType), 'N', num2str(N), 'sample', num2str(sample), '.mat');
            load (fullfile(strcat(loadPath, '\S3C\'),  nameF));
            misS3C(end+1) = missrateS3C;%missrate(end);
        end
        S3C_list(cnt,1) = mean(misS3C);
    end
    
    figure, hold on
    c = 'rgbmkcrgbmkcrgbmkcrgbmkcrgbmkcrgbmkc';
    x = 40:10:90;
    plot (x,SCC_list, 'color', c(1),'marker', '*'),str{1}= 'SSC';
    plot (x,S3C_list, 'color', c(2),'marker', '*'),str{2}= 'S^3C';
    plot (x,probSCC_list, 'color', c(3),'marker', '*'),str{3}= 'Prob SSC';
    grid on
    
    legend(str)
    xlabel('%Intersect','fontsize', 20)
    ylabel('%Error','fontsize', 20)
    xlim([40,90])
    
end
k = 1;