clear 
close all

%% Plotting (20220911) Oscillating force data
load('P:\projects\human\VOP STIM\DATA\Chronic_DBS\BP\20220923\chronic_20210914_baseline_v2\param.mat')
F(1) = figure;

oscTITLE = {'55Hz 35-55%-long', '55 Hz 55-70%-long', '55Hz 55-70%-short','55Hz 35-55%-short','55Hz 35-55%-long', '55 Hz 55-70%-long', '55Hz 55-70%-short','55Hz 35-55%-short','NO STIM 35-55%-long', 'NO STIM 55-70%-long', 'NO STIM 55-70%-short','NO STIM 35-55%-short'};

oscIdx = [2  10]; % High 35%-55%
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


%% Calculate the Power Spectrum for the oscillating task

% close all
trialNames = {'55Hz 55-70%-long','NO STIM 55-70%-long'};
rms = [];
rmsFilt = [];
axes1 = [];
axes2 = [];
axes3 = [];
filtreal = {};

traces = {};
Fs = 60;
window = 500; %in msec


for aa = 1:length(oscIdx)
    trial = filter_data(param.trial(oscIdx(aa)).force);
    idx = 1;
    count_odd = 1;
    count_even = 1;

    for ii = 1:20
%      for ii = 1:length(trial)/(Fs*window/1000)
%         traces{aa,1}(ii,:) = trial(idx:idx+(window/1000*Fs)-1);
% 
%         idx = idx+(window/1000*Fs);
        if rem(ii,2) == 1
            traces{aa,1}(count_odd,:) = trial(idx:idx+(2*Fs)-1);
            idx = idx+(2*Fs);
            count_odd= count_odd+1;
        else
            holdSection = trial(idx:idx+(4*Fs)-1);
            idx_hold = 1;
            for tt = 1:2/(window/1000) 
                traces{aa,2}(count_even,:) = holdSection(idx_hold:idx_hold+(window/1000*Fs)-1);
                count_even=count_even+1;

                idx_hold = idx_hold+(window/1000*Fs);
            end
            idx = idx+(4*Fs);
        end
    end
end

T = 1/Fs;

powerMat = {};

figure;
for n = 1:size(traces,1)
    for nn = 2%1:size(traces,2)
         t = 0:1/Fs:window/1000-1/Fs;
%         if nn==1
%             t = 0:1/Fs:1-1/Fs;
%         else
%             t = 0:1/Fs:2-1/Fs;
%         end

        for pp = 1:size(traces{n,nn},1)
            single = traces{n,nn}(pp,:);
            l = length(single);
            f = (0:l-1)*(Fs/l);
            
%             y = fft(single);
%             power = abs(y).^2/l;

            xTable = timetable(seconds(t'),single');
            
            [pxx,freq] = pspectrum(xTable,'FrequencyLimits',[1 30]);

%             figure;
%             plot(freq,pow2db(pxx))
            powerMat{n}(pp,:) = pxx;

        end
    end
    hold on 
    averagePower = mean(powerMat{n},1);
    plot(freq,pow2db(averagePower))
    xlabel('Frequency (Hz)')
    ylabel('Power Spectrum (dB)')
    title('HOLD 500 msec window')
    legend("stim",'no stim')
end





%% Average Force Curve
% close all
axes = [];
allTrace = {};
filtTrace = {};
allSinTrace = {};
rmsMatFilt = [];
rmsMat = {};

for ii = 1:length(oscIdx)

    idx = floor(length(param.trial(oscIdx(ii)).sinForce(1,:))/5);
    st = 1;
    ed = idx;
    x_time = param.trial(oscIdx(ii)).screenTime(1,st:idx);

    filteredTraceWhole = filter_data(param.trial(oscIdx(ii)).force);
    for n = 1:5

        filtSect = filteredTraceWhole(:,st:ed);
        sect = param.trial(oscIdx(ii)).force(:,st:ed);
        sinSect = param.trial(oscIdx(ii)).sinForce(:,st:ed);
        st = st+idx;
        ed = ed+idx;

        allSinTrace{ii}(n,:) = sinSect;
        allTrace{ii}(n,:) = sect;
        
        filtTrace{ii}(n,:) = filtSect; %+ (mean(sinSect(1:60))-mean(filtSect(1:60)));

        rmsMat{ii}(n,:) = sqrt(mean((sect-sinSect).^2));

        rmsMatFilt{ii}(n,:) = sqrt(mean((filtSect-sinSect).^2));
    end
    figure(5)
    ax = subplot(nTrl,1,ii);
    hold on
    plot(x_time, param.trial(oscIdx(ii)).sinForce(1,idx+1:idx*2),'k--','LineWidth',1.5)
    plot(x_time, param.trial(oscIdx(ii)).sinForce(1,idx+1:idx*2)+[.05 -.05]','k:','LineWidth',1.5)
    plot(x_time, allTrace{ii}(1:5,:),'LineWidth',2)
    title(oscTITLE{oscIdx(ii)})
    
    axes(end+1) = ax;
%     
    figure(6)
    ax1 = subplot(nTrl,1,ii);
    hold on
    plot(x_time, param.trial(oscIdx(ii)).sinForce(1,idx+1:idx*2),'k--','LineWidth',1.5)
    plot(x_time, param.trial(oscIdx(ii)).sinForce(1,idx+1:idx*2)+[.05 -.05]','k:','LineWidth',1.5)
    plot(x_time, filtTrace{ii}(1:5,:),'LineWidth',2)
    title([oscTITLE{oscIdx(ii)} ' Filtered'])
    
    axes1(end+1) = ax1;

end
linkaxes(axes,'xy')
linkaxes(axes1,'xy')


%% Boostrap Significance

[ci95,rejectNull] = bootstrapCompMeans(rmsMatFilt{2}(:,:),rmsMatFilt{1}(:,:),10000,.05)
% [ci95,rejectNull] = bootstrapCompMeans(rmsMatFilt{1}(:,:),rmsMat{2}(:,:),10000,.05)


%% PLOT RMS 

figure;

xErrors = [1 2 3 4];

errorVals = [rmsMat{1}(1:5)'; rmsMat{2}(1:5)'; rmsMatFilt{1}(1:5)';rmsMatFilt{2}(1:5)'];

bar(mean(errorVals'))

stError = std(errorVals')/(5^0.5);

hold on
errorbar(xErrors,mean(errorVals'),stError,'k','linestyle','none')
j = swarmchart(xErrors, errorVals',1000,'red','.');

hold off
trialNames = {'55Hz 55-70%-long','NO STIM 55-70%-long', 'STIM FILT', 'NO STIM FILT'};
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