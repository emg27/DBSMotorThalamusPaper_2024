clear 
close all

%% Plotting (20220911) Oscillating force data
% load('P:\projects\human\VOP STIM\DATA\Chronic_DBS\BP\20220911\chronic_20220911_50Hz\param.mat')%'F:\RNEL\HUMAN\VOP\20220211\subject05_20220211\param_oscillation.mat')


load('P:\projects\human\VOP STIM\DATA\Chronic_DBS\BP\20230630\chronic_20210914_baseline_v2\param.mat')
F(1) = figure;

oscTITLE = {'55Hz 50-70%-short', '55 Hz 35-55%-short', '130Hz 55-70%-short','No Stim 55-70%-short', '55Hz 55-70%-short'};

oscIdx = [5 4]; % High 35%-55%
% oscIdx = [2 5 7]; % Low 20-35%
time = param.trial(oscIdx(1)).screenTime(1,:) - param.trial(oscIdx(1)).screenTime(1);
nTrl = length(oscIdx);
for n = 1:nTrl

    time = param.trial(oscIdx(n)).screenTime(1,:) - param.trial(oscIdx(n)).screenTime(1);
    subplot(nTrl,1,n), hold on
    plot(time, param.trial(oscIdx(n)).sinForce(1,:),'k--','LineWidth',1.5)
    plot(time, param.trial(oscIdx(n)).sinForce(1,:)+[.05 -.05]','k:','LineWidth',1.5)
    plot(time, param.trial(oscIdx(n)).force,'LineWidth',2)
    title(oscTITLE{oscIdx(n)})
    %ylim([-0.05 0.35])
end
xlabel('Time (s)'),ylabel('Normalized Force'),legend('Target Force','Upper Bound','Lower Bound',"Trace")

F(1).Name = '20220911_OscillatingForcePlot';

%% Average Force Curve
% close all
axes = [];
allTrace = {};
filtTrace = {};
allSinTrace = {};
rmsMatFilt = [];
rmsMat = {};


for ii = 1:length(oscIdx)

    idx = floor(length(param.trial(oscIdx(ii)).sinForce(1,:))/3);
    st = 1;
    ed = idx;
    x_time = param.trial(oscIdx(ii)).screenTime(1,st:idx);

%     filteredTraceWhole = filter_data(param.trial(oscIdx(ii)).force);
    count = 1;
    for n = 1:3

%         filtSect = filteredTraceWhole(2,st:ed);
%         sect = mean(param.trial(oscIdx(ii)).force(:,st:ed),1);
%         sinSect = mean(param.trial(oscIdx(ii)).sinForce(:,st:ed),1);

        sect = param.trial(oscIdx(ii)).force(:,st:ed);
        sinSect = param.trial(oscIdx(ii)).sinForce(:,st:ed);
        st = st+idx;
        ed = ed+idx;
        for tt = 1:size(sect,1)
            allSinTrace{ii}(count,:) = sinSect(tt,:);
            allTrace{ii}(count,:) = sect(tt,:);
            
    %         filtTrace{ii}(n,:) = filtSect; %+ (mean(sinSect(1:60))-mean(filtSect(1:60)));
    
            rmsMat{ii}(count,:) = sqrt(mean((sect(tt,:)-sinSect(tt,:))).^2);
            count = count+1;
    
    %         rmsMatFilt{ii}(n,:) = sqrt(mean((filtSect-sinSect).^2));
        end
    end
    figure(5)
    ax = subplot(nTrl,1,ii);
    hold on
    plot(x_time, param.trial(oscIdx(ii)).sinForce(1,idx+1:idx*2),'k--','LineWidth',1.5)
    plot(x_time, param.trial(oscIdx(ii)).sinForce(1,idx+1:idx*2)+[.05 -.05]','k:','LineWidth',1.5)
    plot(x_time, allTrace{ii}(:,:),'LineWidth',2)
    title(oscTITLE{oscIdx(ii)})
    
    axes(end+1) = ax;

end


%% Boostrap Significance

% [ci95,rejectNull] = bootstrapCompMeans(rmsMatFilt{2}(:,:),rmsMatFilt{1}(:,:),10000,.05)
% [ci95,rejectNull] = bootstrapCompMeans(rmsMat{1}(:,:),rmsMat{2}(:,:),10000,.05)
[ci95,rejectNull] = bootstrapCompMeans(rmsMat{1}(:,:),rmsMat{2}(:,:),10000,.05)
% [ci95,rejectNull] = bootstrapCompMeans(rmsMat{2}(:,:),rmsMat{3}(:,:),10000,.05)

%% PLOT RMS 

figure;

xErrors = [1 2];

errorVals = [rmsMat{1}(:)'; rmsMat{2}(:)'];

bar(mean(errorVals'))

stError = std(errorVals')/(5^0.5);

hold on
errorbar(xErrors,mean(errorVals'),stError,'k','linestyle','none')
j = swarmchart(xErrors, errorVals',1000,'red','.');

hold off
trialNames = {oscTITLE{oscIdx}};
set(gca,'xticklabel',trialNames)

%% FUNCTIONS

function filtered_LFP = filter_data(LFP)
fs_chan = 60;

filterbands_Line=[1 4];
[z, p, k] = butter(1, filterbands_Line/(fs_chan/2), 'stop');
[sos_line,g_line] = zp2sos(z, p, k, 'down', 'two');

% filterbands=[15,5000];
% [z,p,k] = butter(2, filterbands/(fs_chan/2), 'bandpass');
% [sos,g] = zp2sos(z, p, k, 'down', 'two');


% filtered_LFP=filtfilt(sos,g,double(LFP'))';
filtered_LFP=filtfilt(sos_line,g_line,double(LFP'))';
end   

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
    ci95 = quantile(diffSampMeans, [alpha/2, 1-(alpha)]);
    
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