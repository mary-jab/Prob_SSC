function [sPath] = plotFigure (QMat,ZKSym, missrate, dataInfoOpt, iteration, errorPre, isSave)
sPath = [];
if nargin<4
    iteration = 0;
end
missrate = round(missrate*100, 2);
errorPre = round(errorPre,2);
if nargin <5
    if isfield(dataInfoOpt, 'errorPrc')
        ttl = ['Intersect: ' num2str(dataInfoOpt.errorPrc) ', MissClassification: ' num2str(missrate)];
    else
        ttl = ['class: ' num2str(dataInfoOpt.k) ', MissClassification: ' num2str(missrate)];
    end
elseif isfield(dataInfoOpt, 'errorPrc')
    ttl = (['Intersect: ' num2str(dataInfoOpt.errorPrc) ', MissClassification: ' num2str(missrate) ', Nan: ' ((isnan(errorPre)))]);
else
    ttl = (['Class: ' num2str(dataInfoOpt.k) ', MissClassification: ' num2str(missrate) ', Nan: ' num2str(errorPre)]);
end
figure;
subplot(121), imagesc(abs(QMat)*800);colorbar;
for i=1:size(ZKSym,1)
zz(i,:) = ZKSym(i,:)/norm(ZKSym(i,:),2);
end
subplot(122), imagesc((zz)*800);colorbar;
set(gcf, 'Position', [400, 400, 1200, 400])

title(ttl);

%% save
if isSave
    if isfield(dataInfoOpt, 'errorPrc')
        sPath = [dataInfoOpt.savePath '/Details/'];%/ErrPrc' num2str(dataInfoOpt.errorPrc)];
        mkdir(sPath);mkdir([sPath '/fig']);
        saveas(gcf,[sPath '/fig/detail_' num2str(dataInfoOpt.k) '_Err_' num2str(dataInfoOpt.errorPrc) ...
            '_smpl_'  num2str(dataInfoOpt.itt) '_itr_' num2str(iteration) '.fig' ])
    else
        sPath = [dataInfoOpt.savePath '/Details/class' num2str(dataInfoOpt.k)];mkdir(sPath);mkdir([sPath '/fig']);
        saveas(gcf,[sPath '/fig/detail_' num2str(dataInfoOpt.k)  ...
            '_sample_'  num2str(iteration) '.fig' ])
    end
end
end




