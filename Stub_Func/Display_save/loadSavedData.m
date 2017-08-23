clear all, close all
for cls = 2:4
    saveArr = []; cnt = 0; S3C_list = [];
    for noise = 40/100: 10/100:90/100
        %         loadPath = strcat(pwd,'\savedRes\GuassianNois\Cls_15\Subspace_noise_', num2str(noise), '\ambiant100');
        %         N = 150;
        cnt = cnt+1;
%         loadPath = strcat(pwd,'/savedRes/Intersect/cls' ,num2str(cls) ,'/Subspace_noise_', num2str(noise), '/ambiant200');
%         N = cls*100;
        cls = 2;
        loadPath = strcat(pwd,'\savedRes\\YaleBCrop025\Class_', num2str(cls));
        N = cls * 64;
        normType = 2;
        
        % misArr = zeros(20,10);
        for sample = 1:160
            nameF =strcat('normType', num2str(normType), 'N', num2str(N));
            nameF =strcat(nameF, 'sample', num2str(sample));
            nameF =strcat(nameF, '.mat');
            load (fullfile(loadPath,  nameF));
            misArr(sample, :) = [missrate, ones(1, 10-length(missrate))* missrate(end)];
            misEnd(sample) = missrate(end);
        end
        
        sz = 70;
        [a,  idx] = sort(misEnd);
         misEnd(idx(sz:end)) = [];
         misArr((idx(sz:end)),:) = [];
        idxList(cnt, :) = idx(1:sz);
        
        probSCC_list(cnt, 1) = mean(misEnd)*100;
        %         SCC_list(cnt,1) = mean(misArr(:,1))*100;
        
        misS3C = []; SSC  = [];
        for sample = idxList(cnt,:)
            nameF =strcat('normType', num2str(normType), 'N', num2str(N), 'sample', num2str(sample), '.mat');
            load (fullfile(strcat(loadPath, '/S3C/'),  nameF));
            misS3C(end+1) = missrateS3C;%missrate(end);
            SSC(end+1) = 1-eval_iter(1);
            a = find(eval_iter>0);
            eval_iter(eval_iter==0) = eval_iter(a(end));
            s3cArrmisArr(sample, :) = 1-eval_iter;
        end
        SCC_list(cnt,1) = mean(SSC)*100;
        S3C_list(cnt,1) = mean(misS3C)*100;
        
        figure, hold on, set(gca,'fontsize',18),
        c = 'krbmkcrgbmkcrgbmkcrgbmkcrgbmkcrgbmkc';
        x = 1:10;
        str = [];
        plot (x,mean(s3cArrmisArr), ':b+', 'LineWidth' , 2),str{1}= 'S^3C';
        plot (x, mean(misArr), '-r*', 'LineWidth' , 2),str{2}= 'Prob SSC';
        grid on
        
        %     xlabel('%Intersect','fontsize', 20)
        %     ylabel('%mis-class','fontsize', 20)
        legend(str)
        %     xlim([40,90])
        %     ylim([0 55])
        %     title(['#Classes: ' num2str(cls)])
        
        
        
        
    end
    %     tmp = S3C_list(end, end); S3C_list(end, end) = SCC_list(end,end); SCC_list(end,end) = tmp;
    figure, hold on, set(gca,'fontsize',18),
    c = 'krbmkcrgbmkcrgbmkcrgbmkcrgbmkcrgbmkc';
    x = 40:10:90;
    str = [];
    plot (x,SCC_list, '--k+', 'LineWidth' , 2),str{1}= 'SSC';
    plot (x,S3C_list, ':b+', 'LineWidth' , 2),str{2}= 'S^3C';
    plot (x,probSCC_list, '-r*', 'LineWidth' , 2),str{3}= 'Prob SSC';
    grid on
    
    xlabel('%Intersect','fontsize', 20)
    ylabel('%mis-class','fontsize', 20)
    legend(str)
    xlim([40,90])
    ylim([0 55])
    title(['#Classes: ' num2str(cls)])
    
end
k = 1;