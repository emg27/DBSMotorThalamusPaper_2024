clear 
close all
%%
fs=30000;

%DEFINE FILTERS
filterbands_Line=[58,62];
[z, p, k] = butter(2, filterbands_Line/(fs/2), 'stop');
[sos_line,g_line] = zp2sos(z, p, k, 'down', 'two');


%% Load

% PATHNAME = 'P:\projects\human\VOP STIM\DATA\Chronic_DBS\TBI01\20220911\';
% trialsToPlot = [3 12 17];
% trialNames = {'55Hz','80Hz','noStim'};

% PATHNAME = 'P:\projects\human\VOP STIM\DATA\Chronic_DBS\TBI01\20220826\';
% trialsToPlot = [5 11];
% trialNames = {'55Hz','noStim'};
% 
% PATHNAME = 'P:\projects\human\VOP STIM\DATA\Chronic_DBS\TBI01\20220511\';
% trialsToPlot = [12 14 17];
% trialNames = {'50Hz','noStim','50Hz'};

% PATHNAME = 'P:\projects\human\VOP STIM\DATA\Chronic_DBS\TBI01\20220914\';
% trialsToPlot = [5 11 20];
% trialNames = {'noStim','noStim','50Hz'};

% PATHNAME = 'P:\projects\human\VOP STIM\DATA\Chronic_DBS\TBI01\20220919\';
% trialsToPlot = [4 12 22];
% trialNames = {'55Hz','55Hz','noStim'};

PATHNAME = 'P:\projects\human\VOP STIM\DATA\Chronic_DBS\TBI01\20220923\';
trialsToPlot = [4 8 14 18 20 24];
trialNames = {'55Hz','55Hz','55Hz','55Hz','noStim','noStim'};


% PATHNAME = 'P:\projects\human\VOP STIM\DATA\Chronic_DBS\TBI01\20221017\';
% trialsToPlot = [2 6 10 14 19 23];
% trialNames = {'55Hz', '55Hz', 'noStim','noStim', '55Hz','55Hz'};

% PATHNAME = 'P:\projects\human\VOP STIM\DATA\Chronic_DBS\TBI01\20230630\';
% trialsToPlot = [6 10 15];
% trialNames = {'55Hz', '130Hz', 'noStim'};


for t=1:length(trialsToPlot)
    
    trialnumber = trialsToPlot(t);
    FILENAME = sprintf('datafile%04d',trialnumber);
    [NS5,~,~] = loadExperimentData(FILENAME,PATHNAME);
    
    analogCh = find(contains({NS5.ElectrodesInfo.Label},'analog 1'));

    disp(['File ' num2str(trialnumber)]);
    data{t} = NS5.Data(analogCh(1),:);

    dataFilt{t}=filtfilt(sos_line,g_line,double(data{t}'))';

end

%% Selecting MVC

manualLabel = false;

close all
maxVals = [];
MVC = [];
axes = [];

maxValCell = {};
maxValMean = [];

for t=1:length(trialsToPlot)
    maxVal = [];
    if manualLabel
            figure;

            plot(dataFilt{t})

            title(trialNames{t})

            [x,y,button]=ginput;
            windowTimes{t} = [x,y,button];
%        
            save([PATHNAME 'MVC_startStopInds.mat'],'windowTimes')
            close all
    else
            windowTimes = load([PATHNAME 'MVC_startStopInds.mat']).windowTimes;
    end
   
    for pks = 1:length(windowTimes{t})
        startT = windowTimes{t}(pks,1);

        if pks+1<length(windowTimes{t})
            endT = windowTimes{t}(pks+1,1);
        else
            endT = length(dataFilt{t});
        end

        window = dataFilt{t}(startT:endT);
        maxVal(end+1) = max(window);
        

    end
    maxValMean(t) = mean(maxVal);
    maxValCell{t} = maxVal;
end

%% Plot and Stats
figure;
maxMat=[];
for ii = 1:length(trialsToPlot)
    maxMat = [maxMat maxValCell{ii}'];
end

bar(mean(maxMat))
xticklabels(trialNames)

xErrors = [1:length(trialsToPlot)];
stError = std(maxMat)/(size(maxMat,1)^0.5);

hold on
errorbar(xErrors,mean(maxMat),stError,'k','linestyle','none')
j = swarmchart(xErrors, maxMat',1000,'red','.');

hold off

 trialMat = [1 2
             ];

        alpha = .05;
        offset = 0;
        for comb = 1:size(trialMat,1)
            trial1 = trialMat(comb,1);
            trial2 = trialMat(comb,2);
    
            [CI,sig]=bootstrapCompMeans(maxValCell{trial1},maxValCell{trial2},10000,alpha/size(trialMat,1));
    
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

end