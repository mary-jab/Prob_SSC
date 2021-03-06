function [Y, s, savePath] = ReadData( input, opt)

if (strcmp(input{1}, 'syn') && strcmp(input{2}, 'GuassianNoise'))
    inputPath = ['../Generate_SyntheticData/genData/GuassianNoise_Cls_' num2str(opt.clas) '/'];
    fileName =  [ inputPath 'Subspace_GuassNoise_Cls_' num2str(opt.clas) '_dim_' ...
        num2str(opt.dim) '_ambient_' num2str(opt.amb) '_noisePrc_' num2str(opt.noise) '_' num2str(opt.sample) '.mat'];
    load (fileName);
    [n, N1, numGrps] = size(Y);
    N= N1*numGrps;
    Y = reshape(( Y(:,:,1:numGrps)),n, N);
    savePath = ['savedRes/GuassianNois/Cls_' num2str(opt.clas) '/Subspace_noise_' num2str(opt.noise) '/ambiant' num2str(opt.amb) '/'];
    
elseif (strcmp(input{1}, 'syn') && strcmp(input{2}, 'intersect'))
    inputPath = ['../Generate_SyntheticData/genData/Intersect'  num2str(opt.clas) 'Cls/'];
    fileName =  [inputPath 'Subspace_intersects_sDim_' num2str(opt.dim) '_subspace' num2str(opt.clas)...
        '_ambiant' num2str(opt.amb) '_prc' num2str(opt.noise) '_iter' num2str(opt.sample) '.mat'];
    load (fileName);
    Y = YMat;
    s = groundTruth;%closeCluster;%
    [n, N1, numGrps] = size(Y);
    N= N1*numGrps;
    %     Y = reshape(( Y(:,:,1:numGrps)),n, N);
    savePath = ['savedRes/Intersect/cls'  num2str(opt.clas) '/Subspace_noise_' num2str(opt.noise) '/ambiant' num2str(opt.amb) '/'];
    
elseif (strcmp(input{1}, 'motion'))
    inputPath =  '..\GRealData\Hopkins155\';
    d = dir(inputPath);
    cnt = 0;
    for i = 1:length(d)
        if ( (d(i).isdir == 1) && ~strcmp(d(i).name,'.') && ~strcmp(d(i).name,'..') )
            cnt = cnt+1;
            if (cnt == opt.sample)
                filepath = d(i).name;
                f = dir([inputPath filepath]);
                foundValidData = false;
                for j = 1:length(f)
                    if ( ~isempty(strfind(f(j).name,'_truth.mat')) )
                        ind = j;
                        foundValidData = true;
                        break
                    end
                end
                eval(['load ' [inputPath filepath '\' f(ind).name]]);
                N = size(x,2);
                F = size(x,3);
                D = 2*F;
                Y = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
                break;
            end
        end
    end
    opt.clas = max(s);
    savePath = ['savedRes/Hopkins/cls'  num2str(opt.clas) '/' filepath '/'];
    
elseif (strcmp(input{1}, 'handwritten'))
    if ~isdeployed
        inputPath = '../GRealData/UCI_MultipleFeaturesDataSet/mfeat/';
    else
        [inputPath,name]=fileparts(which('mfeat-fou'));
        %         inputPath = opt.pathstr;
    end
    load([inputPath '/mfeat-fou']);
    load([inputPath '/mfeat-kar']);
    load([inputPath '/mfeat-fac']);
    %     load('C:\Users\Student\Dropbox\UCF\UCI_MultipleFeaturesDataSet\mfeat\mfeat-pix')
    %     load('C:\Users\Student\Dropbox\UCF\UCI_MultipleFeaturesDataSet\mfeat\mfeat-zer')
    %     load('C:\Users\Student\Dropbox\UCF\UCI_MultipleFeaturesDataSet\mfeat\mfeat-mor')
    Y  = [mfeat_fou mfeat_fac mfeat_kar]';% mfeat_mor mfeat_pix mfeat_zer]';
    numCls = 10;
    classSz = 40;
    totalSz = 200;
    idx = false(1, totalSz*numCls);
    for i=1:numCls
        cur_idx = (i-1)*totalSz+1;
        idx(cur_idx:cur_idx+classSz-1) = 1;
    end
    
    
    Y = Y(:,idx);
%     clsSize = size(Y,2)/numCls;
    s = [];
    for i=1:numCls
        s = [s, ones(1, classSz)*i];
    end
    savePath = ['savedRes/UCI_MultipleFeaturesDataSet/numCls' num2str(numCls) 'clsSize' num2str(classSz) '_featureSize' num2str(size(Y,1)) '/'];
    
end

end


