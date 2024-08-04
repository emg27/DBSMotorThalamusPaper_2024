clear 
close all

%% INPUT VARIABLES

TrialsToLoad=[1 5 6 3]; % pre-lesion trials
% TrialsToLoad=[2 120 121]; %post-lesion trials

pre = -5;
post = 25;

pathName = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-OP M63-18 Opal\20211208 - Experiments\DATA\mat\';

channels = [1:20];
analog = 3;

%% Load Data

for t=1:length(TrialsToLoad)
    ns5_filename = sprintf('NS5_PAD%04d',TrialsToLoad(t));
    load([pathName ns5_filename '.mat']);

%% Grab timestamps
    StimThresh=16000;
    Fs = 30000;    
    
    stim=findAnalogStim(NS5,StimThresh);
    StimTimes = stim;
    vop = stim;
    StimTimes=StimTimes(StimTimes(:,1)==analog,2);
    vop = vop(vop(:,1)==2,2);
    
    if length(unique(stim(:,1))) > 1
        StimTimes = StimTimes(StimTimes>vop(1));
    end
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
    for aa = 1:size(single,1)
        for ii = 1:size(single(aa).indTrl,1)
            mep = single(aa).indTrl(ii,:);
            mep = mep-mean(mep(1:50));

                %calculate peak to peak
                p2p{t}(aa,ii) = peak2peak(mep);
                
                % calculate envelope and AUC
                env = abs(mep);
                area = trapz(time,env);
                auc{t}(aa,ii) = area;
            
        end
    end
end

%% Plot Envelopes and Boxplots
% close all

emgName = {'APL','FDM','EDC','ECR','FCR','BIC','FDC','Jaw',...
        'Sup. Lip','Cheek'};
evnToPlot = [1:10];

% Hz = {'IC alone'};
Hz = {'IC alone', '50 Hz', '80 Hz','100Hz'};
% Hz = {'IC PRE', '5ms','10ms' '50Hz','80Hz', '100 Hz'};
% Hz = {'IC PRE', 'IC POST', '100Hz POST'};

for aa = evnToPlot
    muscle = emgName{aa};
    p2pMat = [];
    p = [];
    aucMat = [];
    a = [];
    figure(aa)
    for t = 1:length(TrialsToLoad)
        singleP2P = p2p{1,t};
        singleAUC = auc{1,t};
        
        p2pCell{t} = (singleP2P(aa,:)');      
        p2pMat = [p2pMat; p2pCell{t}];
        p = [p; ones(size(p2pCell{t}))*t];

        aucCell{t} = (singleAUC(aa,:)');      
        aucMat = [aucMat; aucCell{t}];
        a = [a; ones(size(aucCell{t}))*t];
    end

    subplot(2,2,1)
    boxplot(p2pMat,p,'Labels',Hz,'symbol','+r')
    title(string(muscle)+' P2P')

    subplot(2,2,2)
    boxplot(aucMat,a,'Labels',Hz,'symbol','+r')
    title(string(muscle)+' AUC')    

    % boostrapping for significance
    trialMat = [1 2
                2 3
                ];

    alpha = .05;
    offset = 0;
    for comb = 1:size(trialMat,1)
        trial1 = trialMat(comb,1);
        trial2 = trialMat(comb,2);

        [CI,sig]=bootstrapCompMeans((auc{1,trial1}(aa,:)),(auc{1,trial2}(aa,:)),10000,alpha/size(trialMat,1));

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

% Percent Variation
emgName = {'APL','FDM','EDC','ECR','FCR','BIC','FDC','Jaw',...
        'Sup. Lip','Cheek'};
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
    title('MK-OP ' + string(muscle)+' P2P')
    xticklabels(Hz)
    ylabel("Percent Variation")

    subplot(2,2,4)
    plot((meanauc(:,aa)-baselineauc)/baselineauc*100,'-o')
    title('MK-OP ' +string(muscle)+' AUC')
    xticklabels(Hz)
    ylabel("Percent Variation")

end

%% COMBINED PERCENT INCREASE
%close all
parts = {'finger','hand','arm','face'};
dataSet = {};
percentSet = {};
figure;
for t = 1:length(parts)
    
    if t==1
        ind = [1 2];
    elseif t==2
        ind = [3 4 5 7];
    elseif t==3
        ind = 6;
    elseif t==4
        ind = [8 9 10];
    end
    baselinep2p = mean(meanp2p(1,ind));
    baselineauc = mean(meanauc(1,ind));
    
    dataSet{t} = meanauc(:,ind);
    
    for tt = 1:length(ind)
        baselineauc = dataSet{t}(1,tt);
        averageVals = dataSet{t}(:,tt);

        averageVals(averageVals < baselineauc) = NaN;
        percentSet{t}(:,tt) = (averageVals-baselineauc)/baselineauc*100;

    end

    subplot(2,4,t)
    plot((mean(meanauc(:,ind),2,'omitnan')-baselinep2p)/baselinep2p*100,'-o')
    title(['MK-OP ' parts{t} ' P2P'])
    xticklabels(Hz)
    ylabel("Percent Variation")

    subplot(2,4,t+4)
    plot(mean(percentSet{t},2,'omitnan'),'-o')
    title(['MK-OP ' parts{t} ' AUC'])
    xticklabels(Hz)
    ylabel("Percent Variation")
    
end

%% PLOT EMG
% close all
% 
% CLEAN TRACES
dataClean = data_sta;
for t=1:length(TrialsToLoad)    
    for aa = 1:length(evnToPlot)
        single = p2p{1,t}(aa,:);
        dataClean{1,t}(aa).avg = nanmean(dataClean{1,t}(aa).indTrl);
    end
end

muscles = [1:10];%[2];
sr_ns5 = 30000;
window = (pre*(sr_ns5/1000)):1:(post*(sr_ns5/1000));

figure;

color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[.7 0 0]};

cAvg = [0 0 0
        1 0 0
        0 1 0];

for chn = muscles
    figure;
    axes = [];
    count = 1;
    
    for t = 1:length(TrialsToLoad)
        
        rawEmg = dataClean{1,t}(chn).indTrl;
        allSTA = mean(rawEmg,1);
        x = window/(sr_ns5/1000);
        
        if ~isempty(rawEmg)
            ax = subplot(1,length(TrialsToLoad),count);
            hold on
            plot(window/(sr_ns5/1000),rawEmg(randperm(size(rawEmg,1),40),:),'color',[.7 .7 .7],'LineWidth',.5)
            plot(x,allSTA,'LineWidth',1.5)
            axes(end+1) = ax;
            title([emgName{chn} ' ' Hz{t}])
            xlim([-5 25])
            hold off
        end
        
        count = count+1;
        
    end
%     legend('IC alone', 'VOP 50 Hz', 'VOP 100 Hz')
    linkaxes(axes,'xy')
   
end


%% FUNCTIONS

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

idxEMG = inclChn<=20;

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
    if idxEMG(chn) == 1
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
    data_sta(chn).avg = nanmean(data_sta(chn).indTrl);
end

idx = inclChn<=20;
numEMG = sum(idx)./2;
if numEMG>0
    temp = inclChn(idx==1);
    % EMG Names
    emgName = {'APL','FDM','EDC','ECR','FCR','BIC','FDC','Jaw',...
        'Sup. Lip','Cheek'};
    
    data_sta = processBroadband(data_sta,numEMG,emgName(unique(ceil(temp/2))));
end

%%


[xNum,yNum] = optSubplotLayout(size(data_sta,1));
for chn = 1:size(data_sta,1)
    %subplot(11,14,chn)
    subplot(xNum,yNum,chn),hold on
    plot(window/(sr_ns5/1000),data_sta(chn).avg)

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
    dataClean(2*n-1).indTrl = dataClean(2*n).indTrl- dataClean(2*n-1).indTrl;
    dataClean(2*n-1).avg = nanmean(dataClean(2*n-1).indTrl);
    dataClean(2*n-1).ElectrodesInfo.ConnectorBank = [emgName{n} '_pin'];
end

dataClean(2:2:2*numEMG) = [];
end

function [data]=filter_EMG(data,sr,hp,lp,order)
% 'preprocesso'

% bandwidth filtering-rectification-low pass filtering

%Butterworth 7th order
%Highpass 50 Hz 
cof=hp;
Wn=(cof*2)./sr;
if Wn>1.0
    Wn=0.99;
end
[B,A] = butter(order,Wn,'high'); %Butterworth 7th order
% data_h = fliplr(filter(B,A,fliplr(filter(B,A,data,[],2)),[],2)); % Zero-phase filter
data = fliplr(filter(B,A,fliplr(filter(B,A,data,[],2)),[],2)); % Zero-phase filter     

%Low pass 10Hz
cof=lp;
Wn=(cof*2)./sr;
if Wn>=1.0
    Wn=0.99;
end
[B,A] = butter(order,Wn,'low'); %Butterworth 7th order
% data_hrl = fliplr(filter(B,A,fliplr(filter(B,A,data_h,[],2)),[],2)); % Zero-phase filter
data = fliplr(filter(B,A,fliplr(filter(B,A,data,[],2)),[],2)); % Zero-phase filter
end


%% FUNCTIONS

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

end
