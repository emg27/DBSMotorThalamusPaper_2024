clear
close all


trials = [3 13];
pathname = 'P:\projects\human\VOP STIM\DATA\Chronic_DBS\BP\20230630\';



for a=1:length(trials)
    filename = sprintf('datafile%04d',trials(a));

    data_ns5{a} = loadExperimentData(filename,pathname);
end


%% Extract traces and Frequency Analysis

%Extract and Filter Traces
close all
Fs = 30000;
T = 30; %seconds
labels = {'stim on','stim off'};
trialCh = [1 1];
trialStart = [2351550 3574910 4748910 
              2752730 3938430 5137880]; %20230630


for aa = 1:length(data_ns5)
    for ii = 1:size(trialStart,2)
        wholeTrace{aa,ii} = filter_data(data_ns5{aa}.Data(1,trialStart(aa,ii):trialStart(aa,ii)+T*Fs));
        %     wholeTrace{aa} = (data_ns5{aa}.Data(33,trialStart(aa):trialStart(aa)+T*Fs));
        figure;
        plot(wholeTrace{aa,ii})
        figure(aa*10)
        hold on
        Y = fft(wholeTrace{aa,ii});
        L = length(wholeTrace{aa,ii});
    
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        plot(f,P1)
        hold on
    end
end

% Frequency Analysis
f = {};
averagePower = {};
traces = {};
Fs = 30000;
window = 2000; %in msec
sections = 2; % in seconds

for aa = 1:length(trials)
    for ii = 1:size(trialStart,2)
        trial = wholeTrace{aa,ii};
        idx = 1*Fs;
        count_odd = 1;
        count_even = 1;
        %     traces{aa}(count_even,:) = trial;
    
        while idx+(sections*Fs)-1 < length(trial)
    
            %         traces{aa}(count_even,:) = trial(idx:idx+(sections*Fs)-1);
            idx_hold = 1;
            holdSection = trial(idx:idx+(sections*Fs)-1);
            for tt = 1:sections/(window/1000)
                traces{aa}(count_even,:) = holdSection(idx_hold:idx_hold+(window/1000*Fs)-1);
    
                idx_hold = idx_hold+(window/1000*Fs);
            end
            count_even=count_even+1;
            idx = idx+(4*Fs);
    
        end
    end
end
%% PLOT
close all

fftMat = {};
pspectrumMat = {};
pedoMat = {};
barAverage = [];
barFreq = [];

for n = 1:size(traces,2)
    for ii = 1:size(trialStart,1)
    t = 0:1/Fs:window/1000-1/Fs;

    for pp = 1:size(traces{n},1)
        single = double(traces{n}(pp,:));

        % FFT
        Y = fft(single);
        L = length(single);

        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f{n} = Fs*(0:(L/2))/L;
        fftMat{n}(pp,:) = P1;
    end

    if n == 1
        c= [.7 .7 .7];
        cAvg = [0 0 0];
        shift = 0;
        width = 1;

    else
        c= [0 0 1];
        cAvg = [0 0 .5];
        shift = 1;
        width = .75;
    end

    figure(4)
    hold on
    averagePower{n} = mean(fftMat{n},1);

    %     plot(f{n},fftMat{n}(:,:),'color',c)
    plot(f{n},(averagePower{n}),'color',cAvg,'LineWidth',1.5)

    title("Average Amplitude Spectrum ")
    xlabel("f (Hz)")
    ylabel("Amplitude (N)")
    legend("stim",'no stim')
    xlim([0 12])
    hold off

    skip = 4;
    freqLim = [1 12];

    figure(5)
    hold on

    barFreq(n,:) = f{n}(freqLim(1)*skip+1:skip:freqLim(2)*skip+1);
    barAverage(n,:) = double(averagePower{n}(freqLim(1)*skip+1:skip:freqLim(2)*skip+1));

    b = bar(barFreq(n,:),barAverage(n,:),width,'FaceAlpha',.5,'FaceColor',cAvg);
    hold off

    ax = gca;
    x_ticks = freqLim(1):freqLim(2);

    ax.XTick = x_ticks;
    xlabel("Hz")
    ylabel("Amplitude (N)")
    legend("stim on", "stim off")

    end
end
alpha = .05;
count = 1;
errorStim = [];
errorNoStim = [];

for idx = freqLim(1)*skip+1:skip:freqLim(2)*skip+1
    [CI,sig]=bootstrapCompMeans(fftMat{1}(:,idx),fftMat{2}(:,idx),10000,alpha);

    errorStim(1,count) = std(fftMat{1}(:,idx))/(length(fftMat{1}(:,idx)))^0.5;

    errorNoStim(1,count) = std(fftMat{2}(:,idx))/(length(fftMat{2}(:,idx)))^0.5;

    if sig
        hold on
        plot(x_ticks(count),averagePower{1}(idx)+2,'*k')
        annotation('textbox',[.1 .9 .1 .1],'String',"alpha = "+alpha)
        hold off
    end

    count = count+1;
end

hold on
errorbar(x_ticks,barAverage(1,:),errorStim(1,:),'k','linestyle','none');
errorbar(x_ticks,barAverage(2,:),errorNoStim(1,:),'k','linestyle','none');

%% Boxplot
% close all
statLim = [2 3];
skip = 4;
powerMat = [];
p = [];
powers = {};
seg = {};

for n = 1:length(fftMat)

    window = statLim(1)*skip+1:statLim(2)*skip+1;
%     seg{n} = reshape(fftMat{n}(:,statLim(1)*skip+1:statLim(2)*skip+1),1,[])';
    seg{n} = averagePower{n}(window)';
    powerMat = [powerMat; seg{n}];
    p = [p; ones(size(seg{n}))*n];

end


figure;
boxplot(powerMat,p,'Labels',labels,'symbol','+r')

alpha = .05;
[CI,sig2]=bootstrapCompMeans(seg{1},seg{2},10000,alpha);

if sig2
    yt = get(gca, 'YTick');
    axis([xlim    0  ceil(max(yt)*1.3)])
    annotation('textbox',[.1 .9 .1 .1],'String',"alpha = "+alpha)
    xt = get(gca, 'XTick');
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*1.13, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
    hold off
end
%% FUNCTIONS

function filtered_LFP = filter_data(LFP)
fs_chan = 30000;

% filterbands_Line=[1 6];
% [z, p, k] = butter(1, filterbands_Line/(fs_chan/2), 'stop');
% [sos_line,g_line] = zp2sos(z, p, k, 'down', 'two');

filterbands=[2 5000];
[z,p,k] = butter(1, filterbands/(fs_chan/2), 'bandpass');
[sos,g] = zp2sos(z, p, k, 'down', 'two');


filtered_LFP=filtfilt(sos,g,double(LFP'))';
% filtered_LFP=filtfilt(sos_line,g_line,double(LFP'))';
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

ci95 = quantile(diffSampMeans, [0, 1-alpha]);
%   ci95 = quantile(diffSampMeans, [alpha/2, 1-(alpha/2)]);

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