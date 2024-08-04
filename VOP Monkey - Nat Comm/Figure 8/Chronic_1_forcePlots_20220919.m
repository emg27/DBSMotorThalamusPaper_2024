clear 
close all

%% Plotting (20220911) Oscillating force data
load('P:\projects\human\VOP STIM\DATA\Chronic_DBS\BP\20220919\chronic_20210914_baseline_v2\param.mat')
F(1) = figure;

oscTITLE = {'55Hz 35-55%-short', '55 Hz 55-70%-short', '55Hz 55-70%-long','55Hz 35-55%-long','NO STIM 35-55%-short', 'NO STIM 55-70%-short', 'NO STIM 55-70%-long','NO STIM 35-55%-long'};

% oscIdx = [1 2 5 6]; % short hold
oscIdx = [3 7]; % long hold
% oscIdx = [2 5 7]; % Low 20-35%
time = param.trial(oscIdx(1)).screenTime(1,:) - param.trial(oscIdx(1)).screenTime(1);
nTrl = size(oscIdx,2);
for n = 1:nTrl

    time = param.trial(oscIdx(n)).screenTime(1,:) - param.trial(oscIdx(n)).screenTime(1);
%     subplot(nTrl,1,n), hold on
    plot(time, param.trial(oscIdx(n)).sinForce(1,:),'k--','LineWidth',1.5)
    plot(time, param.trial(oscIdx(n)).sinForce(1,:)+[.05 -.05]','k:','LineWidth',1.5)
    plot(time, param.trial(oscIdx(n)).force,'LineWidth',2)
    title(oscTITLE{oscIdx(n)})
    %ylim([-0.05 0.35])
end
xlabel('Time (s)'),ylabel('Normalized Force'),legend('Target Force','Upper Bound','Lower Bound',"Trace")

F(1).Name = '20220911_OscillatingForcePlot';

% Plotting (20220511) Constant force data
TITLE = {'50 Hz', '80 Hz', 'No Stim'};

%% Calculate the Power Spectrum for the oscillating task

close all
trialNames = {'55Hz 55-70%-long','NO STIM 55-70%-long'};
rms = [];
rmsFilt = [];
axes1 = [];
axes2 = [];
axes3 = [];
filtreal = {};

traces = [];

for aa = 1:length(oscIdx)
    traces(aa,:) = param.trial(oscIdx(aa)).force;
    test = param.trial(oscIdx(aa)).sinForce;
    
    traces(aa+2,:) = filter_data( traces(aa,:));

end

for ii = 1:size(traces,1)
    Fs = 30000;
    T = 1/Fs;
    L = length(traces(ii,:));
    if rem(ii,2) == 1
        t = param.trial(oscIdx(1)).screenTime;  
    else
        t = param.trial(oscIdx(2)).screenTime;
    end
    
%     xTable = timetable(seconds(0:1/60:L)',traces(ii,:)');
    xTable = timetable(seconds(t)',traces(ii,:)');
  
     figure(1)
%      ax1 = subplot(2,1,aa);
    Y = fft(traces(ii,:));

    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    hold on
    f = Fs*(0:(L/2))/L;
    plot(f,P1) 
    title(["Single-Sided Amplitude Spectrum of " trialNames{aa}])
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
    legend('stim on','stim off','filt stim on','filt stim off')
%     axes1(end+1) = ax1;
    hold off

    figure(2)
%     ax2 = subplot(2,1,aa);
    hold on
    [pxx,f] = pspectrum(xTable);
    
    plot(f,pow2db(pxx))
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Power Spectrum (dB)')
    title(['Default Frequency Resolution ' trialNames{aa}])
    legend('stim on','stim off','filt stim on','filt stim off')
%     axes2(end+1) = ax2;

    [P1,f1] = periodogram(traces(ii,:),[],[],Fs,'power');
    hold off
    
   
    figure(3)
     hold on
%     ax3 = subplot(2,1,aa);
    plot(f1,P1)
    grid
    ylabel('P_1')
    title(['Power Spectrum' trialNames{aa}])
    legend('stim on','stim off','filt stim on','filt stim off')
%     axes3(end+1) = ax3;
    hold off


end


[Cxy,f] = mscohere(param.trial(oscIdx(1)).force,param.trial(oscIdx(2)).force,[],[],[],Fs);

thresh = 0.75;
[pks,locs] = findpeaks(Cxy,'MinPeakHeight',thresh);
MatchingFreqs = f(locs);

figure
plot(f,Cxy)
ax = gca;
grid
xlabel('Frequency (Hz)')
title('Coherence Estimate')
ax.XTick = MatchingFreqs;
ax.YTick = thresh;
axis([0 200 0 1])

% linkaxes(axes1,'xy')
% linkaxes(axes2,'xy')
% linkaxes(axes3,'xy')

% figure;
% bar(rms)

%% Average Force Curve
figure;
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
        
        filtTrace{ii}(n,:) = filtSect + (mean(sinSect)-mean(filtSect));

        rmsMat{ii}(n,:) = sqrt(mean((sect-sinSect).^2));

        rmsMatFilt{ii}(n,:) = sqrt(mean((filtSect-sinSect).^2));
    end
    ax = subplot(nTrl,1,ii);
    hold on
    plot(x_time, param.trial(oscIdx(ii)).sinForce(1,idx+1:idx*2),'k--','LineWidth',1.5)
    plot(x_time, param.trial(oscIdx(ii)).sinForce(1,idx+1:idx*2)+[.05 -.05]','k:','LineWidth',1.5)
%     plot(x_time, filtTrace{ii}(1:5,:),'LineWidth',2)
    plot(x_time, allTrace{ii}(1:5,:),'LineWidth',2)
    title(oscTITLE{oscIdx(ii)})
    
    axes(end+1) = ax;
end
linkaxes(axes,'xy')


%% Bootstrap

[ci95,rejectNull] = bootstrapCompMeans(rmsMat{1}(:,:),rmsMat{2}(:,:),10000,.001)

%% Identify the max force
maxF = max(param.trial(oscIdx(1)).sinForce(1,:));
minF = min(param.trial(oscIdx(1)).sinForce(1,:));

% Find the index for each force value
dt = [length(param.trial(oscIdx(1)).sinForce(1,:))-1]/10;
idx = 1:dt:length(param.trial(oscIdx(1)).sinForce(1,:));

dat = repmat(struct('cond',[],'force',nan(2,length(idx)-1),... 
    'rse',nan(3,length(idx)-1)),length(oscIdx),1);

for k = 1:length(oscIdx)
    dat(k).cond = oscTITLE{oscIdx(k)};
    idx = 1:dt:length(param.trial(oscIdx(k)).sinForce(1,:));
    for n = 1:(length(idx)-1)
        tBase = param.trial(oscIdx(k)).sinForce(1,idx(n):idx(n+1)-1);
        tReal = param.trial(oscIdx(k)).force(:,idx(n):idx(n+1)-1);
%         dat(k).rse(:,n) = sqrt(sum((tBase-tReal).^2,2)./3);
        dat(k).rse(:,n) = sqrt(mean((tBase-tReal).^2,2));
        dat(k).force(:,n) = tBase([1 end])';
    end
    for n = 1:2
%         dat(k).combine(:,n) = [dat(k).rse(:,n); dat(k).rse(:,n+4); dat(k).rse(:,n+8);dat(k).rse(:,n+12);dat(k).rse(:,n+16);];
        dat(k).combine(:,n) = [dat(k).rse(:,n); dat(k).rse(:,n+2)];
    end
%     dat(k).combine(1:length(dat(k).rse(:,n)),1) = nan; % Remove the first condition.
end

% Create a boxplot
% TITLES = {'Initial Rise','Late-Stage Rise','Initial Drop','Late-Stage Drop'};
TITLES = {'Rise','Drop'};
F(3) = figure;
for n = 1:2
    subplot(2,2,n)
    boxplot([dat(1).combine(:,n) dat(2).combine(:,n)],{'1','2'})
    title(TITLES{n})
end
F(3).Name = '20220911_OscillatingForceRMSE';

%% PLOT RMS 

figure;

xErrors = [1 2];

errorVals = [rmsMat{1}'; rmsMat{2}'];

bar(mean(errorVals'))

stError = std(errorVals')/(4^0.5);

hold on
errorbar([1 2],mean(errorVals'),stError,'k','linestyle','none')
j = swarmchart(xErrors, errorVals',1000,'red','.');

hold off
trialNames = {'55Hz 55-70%-long','NO STIM 55-70%-long'};
set(gca,'xticklabel',trialNames)



%% FUNCTIONS

function filtered_LFP = filter_data(LFP)
fs_chan = 30000;

% filterbands_Line=[1,15];
% [z, p, k] = butter(2, filterbands_Line/(fs_chan/2), 'stop');
% [sos_line,g_line] = zp2sos(z, p, k, 'down', 'two');

filterbands=[15,5000];
[z,p,k] = butter(2, filterbands/(fs_chan/2), 'bandpass');
[sos,g] = zp2sos(z, p, k, 'down', 'two');


filtered_LFP=filtfilt(sos,g,double(LFP'))';
% filtered_LFP=filtfilt(sos_line,g_line,double(LFP'))';
end   
