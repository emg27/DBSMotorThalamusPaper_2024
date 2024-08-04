clear 
close all

%% INPUT VARIABLES


% TrialsToLoad=[3 12 13 14 15 17 18 19 20 21]; % Delay burst testing
 TrialsToLoad=[175 181 180 182]; % Frequency testing
% TrialsToLoad = [208 210]; % Post lesion data
pre = 5;
post = 25;
analog = 2;


pathName = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-JC M68-22 Johnny Cage\20221005 - Experiment\DATA\mat\';

channels = [129:148];
evnToPlot = [1:10];
emgName = {'EDC', 'ECR', 'FCR', 'FDC', 'ABP', 'FDM', 'Bicep', 'Sup Lip', 'Cheek', 'Jaw'};
% Hz = {'IC alone', 'VOP 50 Hz 400 uS', 'VOP 50 Hz 300 uS', 'VOP 100 Hz'};
Hz = {'IC PRE', '50 Hz', '80Hz','100Hz'};
% Hz = {'IC PRE','IC alone POST', '50 Hz', '80Hz','100Hz'};
% Hz = {'IC pre', '2ms', '5ms', '10ms','50ms','100ms','150ms','300ms','500ms','IC post'};
% Hz = {'IC Pre','50Hz'};


%% Load Data

for t=1:length(TrialsToLoad)
    ns5_filename = sprintf('NS5_PAD%04d',TrialsToLoad(t));
%     NF3_filename = sprintf('NF3_PAD%04d',TrialsToLoad(t));
%     data_ns5 = loadExperimentData(ns5_filename,pathName);
    load([pathName ns5_filename '.mat']);
%     load([pathName NF3_filename '.mat']);

%% Grab timestamps
    StimThresh=16000;
    Fs = 30000;    
    
    stim=findAnalogStim(NS5,StimThresh);
     StimTimes = stim;
    vop = stim;
    StimTimes=StimTimes(StimTimes(:,1)==analog,2);
    vop = vop(vop(:,1)==3,2);
    
%     if length(unique(stim(:,1))) > 1
%         StimTimes = StimTimes(StimTimes>vop(1));
%     end
%% Compute STA
    data_sta{t} = spkTrigAvg(NS5,StimTimes,pre,post,Fs,channels,30000);

end
legend

%% Compute Peak to Peak and AUC
p2p = {};
auc = {};
time = [pre:1/Fs*1000:post];
for t=1:length(TrialsToLoad)
    single = data_sta{1,t};
    for aa = 1:length(single)
        for ii = 1:size(single(aa).indTrl,1)
            mep = single(aa).indTrl(ii,:);
            p2p{t}(aa,ii) = peak2peak(mep);
    
%             env = lowpass(abs(mep),6,30000);
            env = abs(mep);
            area = trapz(time,env);
    
            auc{t}(aa,ii) = area(end);
        end
    end
end

%% Plot Envelopes and Boxplots
close all

for aa = evnToPlot
    muscle = emgName{aa};

    p2pMat = [];
    aucMat = [];
    p = [];
    a = [];

    
    for t = 1:length(TrialsToLoad)
        singleP2P = p2p{1,t};
        singleAUC = auc{1,t};

        p2pCell{t} = (singleP2P(aa,:)');
        aucCell{t} = (singleAUC(aa,:)');
         
        p2pMat = [p2pMat; p2pCell{t}];
        aucMat = [aucMat; aucCell{t}];

        a = [a; ones(size(aucCell{t}))*t];
        p = [p; ones(size(p2pCell{t}))*t];

    end
    
    figure(aa)
    subplot(2,2,1)
    boxplot(p2pMat,p,'Labels',Hz,'symbol','+r','Whisker',1.5)
    title(string(muscle)+' P2P')
    
    subplot(2,2,2)
    boxplot(aucMat,a,'Labels',Hz,'Symbol','+r','Whisker',1.5)
    title(string(muscle)+' AUC')

    trialMat = [1 2];

    alpha = .05;
    offset = 0;
%     for comb = 1:size(trialMat,1)
%         trial1 = trialMat(comb,1);
%         trial2 = trialMat(comb,2);
% 
% %         [CI,sig]=bootstrapCompMeans((p2p{1,trial1}(aa,:)),(p2p{1,trial2}(aa,:)),10000,alpha);
%         [CI,sig]=bootstrapCompMeans((auc{1,trial1}(aa,:)),(auc{1,trial2}(aa,:)),10000,alpha);
% 
%         if sig
%             yt = get(gca, 'YTick');
%             axis([xlim    0  ceil(max(yt)*1.3)])
%             xt = get(gca, 'XTick');
%             hold on
%             plot(xt([trial1 trial2]), [1 1]*max(yt)*(1.1+offset), '-k',  mean(xt([trial1 trial2])), max(yt)*(1.15+offset), '*k')
%             hold off
%         end
%         offset = offset+0.1;
%     end

    

end

% Percent Variation
emgName = {'EDC', 'ECR', 'FCR', 'FDC', 'ABP', 'FDM', 'Bicep', 'Sup Lip', 'Cheek', 'Jaw'};
meanp2p = [];
meanauc = [];

for t = 1:length(p2p)
    for aa = 1:size(p2p{t},1)
        meanp2p(t,aa) = median(p2p{t}(aa,:));
        meanauc(t,aa) = median(auc{t}(aa,:));
    end
end

for aa = 1:size(meanp2p,2)
    baselinep2p = meanp2p(1,aa);
    baselineauc = meanauc(1,aa);
    muscle = emgName{aa};

    figure(aa)
    subplot(2,2,3)
    plot((meanp2p(:,aa)-baselinep2p)/baselinep2p*100,'-o')
    title(string(muscle)+' P2P')
    xticklabels(Hz)
    ylabel("Percent Variation")

    subplot(2,2,4)
    plot((meanauc(:,aa)-baselineauc)/baselineauc*100,'-o')
    title(string(muscle)+' AUC')
    xticklabels(Hz)
    ylabel("Percent Variation")

end

%% COMBINED PERCENT INCREASE
%close all

parts = {'finger','hand','arm'};%,'face'};
figure;
for t = 1:length(parts)
    if t==1
        ind = [5 6];
    elseif t==2
        ind = [1 2 4];
    elseif t==3
        ind = 7;
    elseif t==4
        ind = [];
    end
    baselinep2p = mean(meanp2p(1,ind));
    baselineauc = mean(meanauc(1,ind));
    
    dataSet{t} = meanauc(:,ind);
    for tt = 1:length(ind)
        baselineauc = dataSet{t}(1,tt);
        averageVals = dataSet{t}(:,tt);

        averageVals(averageVals < baselineauc*0.95) = NaN;
        percentSet{t}(:,tt) = (averageVals-baselineauc)/baselineauc*100;

    end

    subplot(2,4,t)
    plot(mean((meanp2p(:,ind)-baselinep2p),2)/baselinep2p*100,'-o')
    title(['MK-JC ' parts{t} ' P2P'])
    xticklabels(Hz)
    ylabel("Percent Variation")

    subplot(2,4,t+4)
    plot(mean(percentSet{t},2,'omitnan'),'-o')
    title([parts{t} ' AUC'])
    xticklabels(Hz)
    ylabel("Percent Variation")
    
end

%% PLOT EMG
% close all
% CLEAN TRACES
dataClean = data_sta;

% for t=1:length(TrialsToLoad)    
%     for aa = evnToPlot
% %         mean_signal=mean(dataClean{1,t}(aa).indTrl,2);
% %         std_signal=std(mean_signal);
% %         idxToRemove=find(abs(mean_signal)>std_signal*2);
% %         dataClean{1,t}(aa).indTrl(idxToRemove,:)=[];
% 
%         single = p2p{1,t}(aa,:);
%         [clean, idxToRemove] = rmoutliers(single);
% %          std_p2p = std(single);
% % 
% %          idxToRemove = find(abs(single)>std_p2p*5);
%         
%         dataClean{1,t}(aa).indTrl(idxToRemove,:) = [];
%         dataClean{1,t}(aa).avg = nanmean(dataClean{1,t}(aa).indTrl,1);
%     end
% end

sr_ns5 = 30000;
window = (pre*(sr_ns5/1000)):1:(post*(sr_ns5/1000));

figure;

color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[.7 0 0]};

cAvg = [0 0 0
        1 0 0
        0 1 0];

for chn = evnToPlot
    figure;
    axes = [];
    count = 1;
    for t = 1:length(TrialsToLoad)

        ax = subplot(1,length(TrialsToLoad),count);
        hold on
        
        title([emgName{chn} ' ' Hz{t}])
        staEmg = dataClean{1,t}(chn).avg;
        rawEmg = dataClean{1,t}(chn).indTrl;
        x = window/(sr_ns5/1000);

%         stdError = std(rawEmg,0,1)/(size(rawEmg,1)^0.5);
        if ~isempty(rawEmg)
        
%             plot(window/(sr_ns5/1000),rawEmg(randperm(size(rawEmg,1),25),:),'color',[.7 .7 .7],'LineWidth',.5)
            plot(window/(sr_ns5/1000),rawEmg(:,:),'color',[.7 .7 .7],'LineWidth',.5)
            plot(x,staEmg,'LineWidth',1.5)
        end

%         curve1 = staEmg + stdError;
%         curve2 = staEmg - stdError;
% 
%         x2 = [x, fliplr(x)];
%         inBetween = [curve1 fliplr(curve2)];
%         fill(x2, inBetween, color{t}, 'faceAlpha',.3)

        axes(end+1) = ax;
        count = count+1;
        xlim([5 25])
    end
%     legend('IC alone', 'VOP 50 Hz', 'VOP 100 Hz')
    linkaxes(axes,'xy')
    hold off
    saveas(gcf,['MK-JC_postLesion_EMGtrace_' emgName{chn}])
end


%% Functions
function [stim] = findAnalogStim(dat_ns5,stimThresh)
% Determine when the stimulations occur. Good Stim Channels: 2, 3, 4.
% Find all the first threshold crossings.
chName = 'analog';
chnIdx = find(contains({dat_ns5.ElectrodesInfo.Label},chName));
inclChn = 1:4;

stimRaw = double(dat_ns5.Data(chnIdx(inclChn),:));
diffStim = diff(stimRaw,[],2);
[stimChn,stimIdx] = find([diffStim]>stimThresh);
stim = [stimChn stimIdx]; % 1st column stimulation channel, 2nd index
stim = sortrows(stim,2);
end

function [data_sta] = spkTrigAvg(dat_ns5,stim,preStim,postStim,sr_ns5,inclChn,F,sr_emg)
% Define the index for each of the stimulation spikes.
if nargin <8
    sr_emg = sr_ns5;
end

window = (preStim*(sr_ns5/1000)):1:(postStim*(sr_ns5/1000));
baseTrue  = 1;

idxEMG = inclChn>=128;

if isempty(inclChn)
    data_sta = repmat(struct('indTrl',nan(size(stim,1),length(window)),...
        'avg',[]),size(dat_ns5.Data,1),1);
    inclChn = 1:size(dat_ns5.Data,1);
else
   data_sta = repmat(struct('indTrl',nan(size(stim,1),length(window)),...
    'avg',[]),length(inclChn),1);
end

figure(1)
for chn = 1:length(inclChn)
    if idxEMG(chn) == 1;
       dat_ns5.Data(inclChn(chn),:) = EMG.filter_EMG(double(dat_ns5.Data(inclChn(chn),:)),sr_emg,50,500,5);
    end
    for n = 1:size(stim,1)
        idx = window + stim(n);
        pos = find(idx>0);
        if idx(end)<= length(dat_ns5.Data)
            if baseTrue
                baseline = mean(dat_ns5.Data(inclChn(chn),stim(n)+[-preStim:-preStim+100]));
                data_sta(chn).indTrl(n,pos) = dat_ns5.Data(inclChn(chn),idx(pos))-baseline;
            else
                data_sta(chn).indTrl(n,pos) = dat_ns5.Data(inclChn(chn),idx(pos));
            end
        else
            temp = dat_ns5.Data(inclChn(chn),idx(pos(1)):end);
            if baseTrue
                baseline = mean(dat_ns5.Data(inclChn(chn),stim(n)+[-preStim:-preStim+100]));
                temp = temp-baseline;
            end
            data_sta(chn).indTrl(n,pos(1)+[1:length(temp)]) = temp;
        end
        data_sta(chn).ElectrodesInfo = dat_ns5.ElectrodesInfo(inclChn(chn));
    end

    data_sta(chn).avg = nanmean((data_sta(chn).indTrl(:,:)));
end

idx = inclChn>=128;
numEMG = sum(idx)./2;
if numEMG>0
    temp = inclChn(idx==1);
    % EMG Names
    emgName = {'EDC', 'ECR', 'FCR', 'FDC', 'ABP', 'FDM', 'Bicep', 'Sup Lip', 'Cheek', 'Jaw'};
    
    data_sta = processBroadband(data_sta,numEMG,emgName(unique(ceil((temp-128)/2))));
end

[xNum,yNum] = optSubplotLayout(size(data_sta,1));
for chn = 1:size(data_sta,1)
    %subplot(11,14,chn)
    subplot(xNum,yNum,chn),hold on
    plot(window/(sr_ns5/1000),data_sta(chn).avg)
%     plot(window,data_sta(chn).avg)
    title([data_sta(chn).ElectrodesInfo.ConnectorBank ...
        num2str(data_sta(chn).ElectrodesInfo.ConnectorPin)])
    xlim([preStim postStim])
end
end

function [dataClean] = processBroadband(data,numEMG,emgName,varargin)
% Clean up the broadband data.\
%numEMG = 10;
sr_ns5 = 30000;
dataClean = data;

% Differential and processing of the EMGs
for n = 1:numEMG
    %temp1 = EMG.filter_EMG(dataClean(2*n).indTrl,sr_ns5,30,500,3);
    %temp2 = EMG.filter_EMG(dataClean(2*n-1).indTrl,sr_ns5,30,500,3);
    %dataClean(2*n-1).indTrl = temp1-temp2;%dataClean(2*n).indTrl- dataClean(2*n-1).indTrl;
    dataClean(2*n-1).indTrl = dataClean(2*n).indTrl- dataClean(2*n-1).indTrl;
    dataClean(2*n-1).avg = nanmean(dataClean(2*n-1).indTrl);
    dataClean(2*n-1).ElectrodesInfo.ConnectorBank = [emgName{n} '_pin'];
end

dataClean(2:2:2*numEMG) = [];
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