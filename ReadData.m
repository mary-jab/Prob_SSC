function [Y, s, savePath] = ReadData( input, opt)

if (strcmp(input{1}, 'syn') && strcmp(input{2}, 'GuassianNoise'))
    inputPath = '../Generate_SyntheticData/genData/GuassianNoise/';
    fileName =  [ inputPath 'Subspace_GuassNoise_Cls_' num2str(opt.clas) '_dim_' ...
        num2str(opt.dim) '_ambient_' num2str(opt.amb) '_noisePrc_' num2str(opt.noise) '_' num2str(opt.samNum) '.mat'];
    load (fileName);
    [n, N1, numGrps] = size(Y);
    N= N1*numGrps;
    Y = reshape(( Y(:,:,1:numGrps)),n, N);
    savePath = ['savedRes/GuassianNoisV/Subspace_noise_' num2str(opt.noise) '/ambiant' num2str(opt.amb) '/'];
    
elseif (strcmp(input{1}, 'syn') && strcmp(input{2}, 'intersect'))
    inputPath = ['../Generate_SyntheticData/genData/Intersect'  num2str(opt.clas) 'Cls/'];
    fileName =  [inputPath 'Subspace_intersects_sDim_' num2str(opt.dim) '_subspace' num2str(opt.clas)...
        '_ambiant' num2str(opt.amb) '_prc' num2str(opt.noise) '_iter' num2str(opt.samNum) '.mat'];
    load (fileName);
    Y = YMat;
    s = groundTruth;
    [n, N1, numGrps] = size(Y);
    N= N1*numGrps;
    %     Y = reshape(( Y(:,:,1:numGrps)),n, N);
    savePath = ['savedRes/Intersect'  num2str(opt.clas) 'Cls/Subspace_noise_' num2str(opt.noise) '/ambiant' num2str(opt.amb) '/'];
    
    
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
    Y = Y(:,1:200*numCls);
    clsSize = size(Y,2)/numCls;
    s = [];
    for i=1:numCls
        s = [s, ones(1, clsSize)*i];
    end
    savePath = ['savedRes/UCI_MultipleFeaturesDataSet/numCls' num2str(numCls) '_featureSize' num2str(size(Y,1)) '/'];
    
end

end


