%% Cortical cord spiking data: raster PSTH heatmap

close all 
clear all
addpath('D:\Documents\ImgProc_Scripts')

%indeces of spinal lead
NeuralData_idx=[33:128];
NChan=length(NeuralData_idx);

%define trials to load
TrialsToLoad=[13];
TrialsHz = {'M1 1 Hz PAD 13','S1 1 Hz','M1 10 Hz PAD 14','S1 10 Hz'};
pathName = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-SC M137-21 Scorpion\20210714 - Experiments\RIPPLE\';
m1Array = [64:-1:33 96:-2:66]-32;
s1Array = [95:-2:65 128:-1:97]-32;

SingleChannelSpikes=[];

PSTH=[];

for t=1:length(TrialsToLoad)
    %% grabs timestamps from NEV
    nev_filename = sprintf('PAD%04d',TrialsToLoad(t));
    load([pathName nev_filename '.mat']);
    nev = NEV;

    electrodes = unique(nev.Data.Spikes.Electrode);
    stim_elec = 5377;%electrodes(end);%
    ind = nev.Data.Spikes.Electrode == stim_elec;
    rawTimes = nev.Data.Spikes.TimeStamp(ind);
    StimTimes = rawTimes(1:3:end);
    
    %% create neural raster plot
    Fs=30000;
    time=[1:nev.MetaTags.DataDuration];
    Raster{t}=zeros(length(NeuralData_idx),length(time));
    
    
    for ii=1:length(nev.Data.Spikes.TimeStamp)
        if find(NeuralData_idx==double(nev.Data.Spikes.Electrode(ii)))
            Raster{t}(nev.Data.Spikes.Electrode(ii)-NeuralData_idx(1)+1, nev.Data.Spikes.TimeStamp(ii))=1;
        end
    end
    
    UtahMapping=[33:128]-32;    
    Raster{t}=Raster{t}(UtahMapping,:);
   
    
    
    %% Computing the PSTH for each channel
    post_stim=60e-3; %window in ms
    pre_stim=15e-3; %pre stim window
    base=10e-3*Fs; %window for baseline in points
    bin=2e-3*Fs; %PSTH bin size in points
    
    time_window=[-pre_stim*Fs:1:post_stim*Fs-1];
    

    for ii=1:length(StimTimes)
        if (StimTimes(ii)+post_stim*Fs-1<length(time))
            for ch=1:NChan            
                SingleChannelSpikes{t,ch}(ii,:)=Raster{t}(ch,[StimTimes(ii)-pre_stim*Fs:StimTimes(ii)+post_stim*Fs-1]);
            %blanking
            SingleChannelSpikes{t,ch}(ii,(pre_stim-3e-3)*Fs:(pre_stim+14e-3)*Fs)=0;
            end
        end
    end
    
    for ch=1:NChan 
    %computing PSTH
        SpikesSum=sum(SingleChannelSpikes{t,ch});
        jj=1;
        for ii=1:floor(length(SpikesSum)/bin)
            PSTH_spikes{t,ch}(ii)=sum(SpikesSum(jj:jj+bin-1)); %spikes count per bin
            %convert to FR: spikes/stims/bin time = spikes/s
            PSTH{t,ch}(ii)=PSTH_spikes{t,ch}(ii);%/(length(StimTimes)/(bin*30000));
            jj=jj+bin;
        end

    %rough baseline correct
        base_bins=floor(base/bin);
        baseline{t,ch}=mean(PSTH{t,ch}(1:base_bins),'all');
        PSTH_base{t,ch}=PSTH{t,ch}-baseline{t,ch};

        spikeBaseline{t,ch} = mean(SingleChannelSpikes{t,ch}(1:base),'all');
        
        spikes_base{t,ch} = SingleChannelSpikes{t,ch} - spikeBaseline{t,ch};

        
    %blanking
        pre_bin = floor((pre_stim-1e-3)*Fs/bin);
        post_bin = floor((pre_stim+14e-3)*Fs/bin);
        PSTH_blank{t,ch}=PSTH_base{t,ch};
        PSTH_blank{t,ch}(pre_bin:post_bin)=0;
    end
    
     
end


%% PLOT
close all
Ntrials=length(TrialsToLoad);

channelsToPlot=[1 10 80 96];
NPlots=length(channelsToPlot);

% Rasters 1:96
figure
for t=1:Ntrials
    for ch=1:NPlots
        subplot(Ntrials,NPlots,ch+(NPlots*(t-1)));
    end
    linkaxes
end
% 
% % PST Histogram
figure
for t=1:Ntrials
    for ch=1:NPlots
        subplot(Ntrials,NPlots,ch+(NPlots*(t-1)))
        time_bins=([0:length(PSTH_base{t,channelsToPlot(ch)})-1]*bin/Fs-pre_stim)*1000;
        bar(time_bins, PSTH_spikes{t,channelsToPlot(ch)},'histc')
        ax=gca;
        ax.YLim=[0 80];
        title(['Trial ' num2str(TrialsToLoad(t)) ' Channel ' num2str(channelsToPlot(ch))]) 
        xlabel('Time (ms)'); ylabel('Firing Rate (Hz)')
        hold on
        plot([0 0],ax.YLim,'r-')
    end
    linkaxes
end

%% Heatmap by channel order
% create custom colormap
custom = zeros(170,3);
left_color = [1 1 1]; % white
mid_color = [1 0 0]; % color according to the attached image
right_color = [1 1 0]; % yellow
redcol = interp1([0, 1], [left_color; mid_color], linspace(0, 1, 120));
%yellcol = interp1([0, 1], [mid_color; right_color], linspace(0, 1, 100));

graycol = linspace(0.5,1,50);

custom(1:50,1) = graycol;
custom(1:50,2) = graycol;
custom(1:50,3) = graycol;
custom(51:170,:) = redcol;%; yellcol];

% plot heatmap
figure
cdata=zeros(48,length(time_bins),Ntrials);%zeros(32,length(time_bins),Ntrials);
arrayCh = [1:48];
elecCount = {};

for t=1:(Ntrials*2)
    subplot(Ntrials*2,1,t)
%     time_bins=([0:length(PSTH{t,1})-1]*bin/Fs-pre_stim)*1000;
    for i=arrayCh
        trialNum = ceil(t/2);
        if mod(t,2)==1
            ii = m1Array(i);
        else
            ii = s1Array(i);
        end
         
        cdata(i,:,t)=PSTH_blank{trialNum,ii}(1:length(time_bins));
        allchSpike{t,i} = cdata(i,:,t);
        elecCount{t}(i,:) = PSTH_blank{trialNum,ii}(1:length(time_bins));
    end
    xvalues = num2cell(time_bins);
    yvalues = num2cell(arrayCh);
    %h = heatmap(xvalues,yvalues,cdata(:,:,t));
    
    
    imagesc(cdata(:,:,t),[-20 50]);
    meanSpike{t} = mean(cdata(:,:,t),1);
    stdSpike{t} = std(cdata(:,:,t))/(length(cdata(:,:,t)))^(.5);
    colormap(custom); 
    colorbar; colormapeditor;
    xticklabels = -10:10:60;
    xticks = linspace(1, length(time_bins), numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    yticklabels = flip([1,16,32,48]);
    yticks = linspace(1, 48, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
    title(['Mk Scorpion - ' TrialsHz{t}]);% Trial ' num2str(TrialsToLoad(t))])
    xlabel('Time (ms)');ylabel('Channel')
%     make_pretty()
end

%% Boostrap
bool = [];
for ii = 1:size(elecCount{1},2)
    [CI,sig] = bootstrapCompMeans(elecCount{1}(:,ii),elecCount{2}(:,ii),10000,.01);
    bool(ii) = sig;
end
%% Average Firing Rate
figure;
for ii = 1:size(meanSpike,2)
    hold on
    x = -11:2:61;
    if ii ==1
        c= [.7 .7 .7];
    else
        c = [.4 .4 .4];
    end
    for aa = 1:size(allchSpike,2)
        data = allchSpike{ii,aa};
        plot(x,data,'color',c)
    end
    color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
    data = meanSpike{ii};
    
    plot(x,data)
    plot(x,bool*30,'*','MarkerSize',10)
%     x = 1:length(data);
    std_dev = stdSpike{ii};
    curve1 = data + std_dev;
    curve2 = data - std_dev;

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