clear
close all

pathName = 'D:\Figures\VOP_NatCom_2024\Figure 8\DLC\';
pathNameInds = 'D:\Figures\VOP_NatCom_2024\Figure 8\';
myFiles = dir(fullfile(pathName, '*.mat'));

trialsToLoad = [1:6];

bodypartlst = {'weight'};
labels = {myFiles.name};
riseTraces = {};

%% Load Data

for t = 1:length(trialsToLoad)
    dlc_filename = myFiles(t).name;
    load([pathName dlc_filename]);
    allData{t} = dat;
    datacol = 2;
    col = 1;
    for ii = 1:length(bodypartlst)

        trial = allData{1,t};
        traces{t}(:,col:col+1) = trial(:,datacol:datacol+1);
        datacol = datacol+3;
        col = col+2;
    end

end

%% Label Start and Peak times

manualLabel = false;

close all
maxVals = [];
MVC = [];
axes = [];

velCell  = {};
smoothCell = {};
p2pCell = {};
maxValCell = {};
maxValMean = [];
colors = {'blue','yellow'};
weight = {'3lbs','8lbs','no weight'};
count=0;

for t=1:length(trialsToLoad)
    %     ax = subplot(1,length(trialsToPlot),t);
    %     axes(end+1)=ax;
    maxVal = [];
    p2pVal = [];
    velVal = [];
    smoothVal = [];
    if manualLabel
%         clear windowTimes
        figure;

        plot(traces{t}(:,2)*-1)
        title(labels{t})

        [x,y,button]=ginput;
        windowTimes{t} = [x,y,button];
        save([pathNameInds 'BP_upperLimb_startStopInds_v2.mat'],'windowTimes')
        close all
    else
        windowTimes = load([pathNameInds 'BP_upperLimb_startStopInds_v2.mat']).windowTimes;
    end
    count = count+1;
    if t==1
        figure;
       
    elseif t==3
        figure;
        count = 1;

    elseif t==5
        figure;
        count = 1;
          
    end
    
    
    for pks = 1:2:length(windowTimes{t})
        startT = floor(windowTimes{t}(pks,1));

        if pks+1 < length(windowTimes{t})+1
            peakT = floor(windowTimes{t}(pks+1,1));
        end

        if pks+1<length(windowTimes{t})
            endT = floor(windowTimes{t}(pks+2,1));
        else
            endT = length(traces{t});
        end

        windowY = traces{t}(startT:endT,2);
       
        if t==1 && ismember(pks,[5 7 9])
            baselineStart = floor(windowTimes{t}(3,1));
%             endT = floor(windowTimes{t}(11,1));
         
            windowY = windowY-traces{t}(baselineStart);
        else
            windowY = windowY-windowY(1);
        end
        windowX = traces{t}(startT:endT,1);
        windowY = -windowY;

        [maxHeight, maxInd] = max(windowY);
        maxVal(end+1) = maxHeight;
        p2pVal(end+1) = peak2peak(windowY);
        velVal(end+1) = windowY(peakT-startT)/(peakT-startT);
%         velVal(end+1) = endT-startT;
%         velVal(end+1) = max(gradient(windowY));

        % Quantify smoothness
        dx = gradient(windowX); % First derivative of x (assuming parameterization along x)
        dy = gradient(windowY); % First derivative of y
        d2x = gradient(dx); % Second derivative of x
        d2y = gradient(dy); % Second derivative of y
        
        % Compute curvature
        curvature = abs(dx .* d2y - dy .* d2x) ./ (dx.^2 + dy.^2).^1.5;

        smoothVal(end+1) = mean(curvature);
        
        riseWindowY = -(traces{t}(startT:peakT,2)-traces{t}(startT));
        riseWindowX = -(traces{t}(startT:peakT,1)-traces{t}(startT));

        hold on
        subplot(1,2,1)
        plot(windowX,windowY,colors{count})
        hold off
        
        hold on
        subplot(1,2,2)
        plot(riseWindowX,riseWindowY,colors{count})
    end

    if t==1
        [val,ind] = min(maxVal);
        maxVal(ind) = [];
        velVal(ind) = [];
    end
  
    
    maxValMean(t) = mean(maxVal);
    maxValCell{t} = maxVal;
    p2pCell{t} = p2pVal;
    velCell{t} = velVal;
    smoothCell{t} = smoothVal;
end

%% Panel H: Plot and Stats

maxMat=[];
maxError = [];

p2pMat = [];
p2pError = [];

weightLab = {'3 pounds', '8 pounds', 'unweighted'};
measurementLab = {'Max Height','P2P','Velocity','Smoothness'};
measurment = 1;

if measurment == 1
    measurementCell = maxValCell;
elseif measurment == 2
    measurementCell = p2pCell;
elseif measurment == 3
    measurementCell = velCell;
elseif measurment == 4
    measurementCell = smoothCell;
end

for ii = 1:length(trialsToLoad)
        maxMat = [maxMat mean((measurementCell{ii}))];
        maxError(ii) = std(measurementCell{ii})/(length(measurementCell{ii})^0.5);
end
figure;
count=0;
axes= [];
for tt = 1:2:length(trialsToLoad)
    count=count+1;
    ax = subplot(1,3,count);
    bar(maxMat(tt:tt+1))
    xticklabels(labels(tt:tt+1))
    title(measurementLab{measurment})
    axes(end+1) = ax;

    xErrors = [1:2];

    hold on
    errorbar(xErrors,maxMat(tt:tt+1),maxError(tt:tt+1),'k','linestyle','none')
    for aa = 1:2
        plot(aa-.1+(rand([1 length(measurementCell{tt+aa-1})])/4),measurementCell{tt+aa-1},'ok')
    end
    hold off

    trialMat = [
        1 2
        ];

    alpha = .05;
    offset = 0;
    for comb = 1:size(trialMat,1)
        trial1 = trialMat(comb,1);
        trial2 = trialMat(comb,2);

        [CI,sig]=bootstrapCompMeans(measurementCell{trial1+(tt-1)},measurementCell{trial2+tt-1},10000,alpha/size(trialMat,1));
        %         [CI,sig]=bootstrapCompMeans(ampMat{trial1},ampMat{trial2},10000,alpha/size(trialMat,1));

        if sig
            yt = get(gca, 'YTick');
            axis([xlim    0  ceil(max(yt)*1.3)])
            annotation('textbox',[.1 .9 .1 .1],'String',"alpha = "+alpha)
            xt = get(gca, 'XTick');
            hold on
                plot(xt([trial1 trial2]), [1 1]*max(yt)*(1.1+offset), '-k',  mean(xt([trial1 trial2])), max(yt)*(1.15+offset), '*k')
            hold off
        end
        offset = offset+0.001;
    end
end
linkaxes(axes,'xy')

% Plot Percent Var
axes1 = [];
values= [];
percentPts = {};
count=0;

figure;
for tt=1:2:length(trialsToLoad)
    count=count+1;

    dataPts = measurementCell{tt};
    baselineMax = maxMat(tt+1);
    stimOnMax = maxMat(tt);
    values = ([values (stimOnMax-baselineMax)/baselineMax*100]);
 
    percentPts{count} = (dataPts-baselineMax)/baselineMax*100;

end
bar(values)
title(measurementLab{measurment})
xticklabels(weightLab)
hold on
for ii = 1:length(percentPts)

    ste = std(percentPts{ii});
    plot(ii,percentPts{ii},'.r',MarkerSize=40)
    errorbar(ii,mean(percentPts{ii}),ste,'k')


end


%% Bootstrapping

function [ci95, rejectNull] = bootstrapCompMeans(dataSet1, dataSet2, bootstrapReps,alpha)

sampMeans1 = nan(1,bootstrapReps);
sampMeans2 = nan(1,bootstrapReps);
diffSampMeans = nan(1,bootstrapReps);
for i=1:bootstrapReps
    % Resample from each dataset with replacement
    bootstrapSamp1 = randsample(dataSet1, length(dataSet1), true);
    bootstrapSamp2 = randsample(dataSet2, length(dataSet2), true);

    % Get means of both samples
    meanSamp1 = mean(bootstrapSamp1);
    sampMeans1(i) = meanSamp1;
    meanSamp2 = mean(bootstrapSamp2);
    sampMeans2(i) = meanSamp2;

    % Get the difference of the means
    diffMeans = meanSamp1 - meanSamp2;
    diffSampMeans(i) = diffMeans;
end

% Calculate confidence interval of difference of means
ci95 = quantile(diffSampMeans, [alpha/2, 1-(alpha/2)]);

% If ci95 contains 0 then don't reject null
if (ci95(1) <= 0) && (ci95(2) >= 0)
    rejectNull = false;
else
    rejectNull = true;
end

%     Plot histogram to confirm
%     figure;
%     hist(diffSampMeans, 100)

end