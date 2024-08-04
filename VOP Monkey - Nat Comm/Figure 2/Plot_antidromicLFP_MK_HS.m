clear all
close all

s1channels = [1:32 33:2:63];
m1channels = [34:2:64 65:96];

%% INPUT VARIABLES
TrialsToLoad=[15 3];
channelsToLoad = [1:96];
brainArea = 'm1'; 
analog = 2; % analog ch 2 or 3 (according to AM stimulator)
pre = 20;
post = 100;

savefile = false;
pathName = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-HS M219-21 Hosu\20220427 - Experiment\DATA\mat\';


%% grabs timestamps from NS5 and LFP from 
StimThresh=16000;
Fs = 30000;

for t=1:length(TrialsToLoad)
    ns5_filename = sprintf('NS5_PAD%04d',TrialsToLoad(t));
    load([pathName ns5_filename '.mat']);
   stimAll=findAnalogStim(NS5,StimThresh);

   StimTimes=stimAll(stimAll(:,1)==analog,2);

% Grabs LFP from all the channels and applies notch filter at 60Hz and
% bandpass at 10-5000Hz
    for iC = 1:length(channelsToLoad) 
        ch = channelsToLoad(iC);
        rawLFP{t,iC} = double(NS5.Data(ch,:));
        filteredLFP{t,iC} = notchFilter_LFP(rawLFP{t,iC});
    end    

%% Build STA
    for iC = 1:length(channelsToLoad)
        singleLFP = filteredLFP{t,iC};
        for aa = 1:length(StimTimes)
           stim = StimTimes(aa);
           if stim+(post*Fs/1000) < length(singleLFP)
               
               lowerbound = stim-(pre*Fs/1000);
               upperbound = stim+(post*Fs/1000);
               window = singleLFP(lowerbound:upperbound);
               %baseline correct
               window = window - mean(window(1:(stim-lowerbound-100))); 
               window(1:667) = 0; % 1ms blanking of stim artifact, IC antidromic

               summed_LFP_single{t,iC}(aa,:) = window;
           end
        end 
    end   

    for iC = 1:length(channelsToLoad)
    %Clean Variables
        summed_LFP_cleaned{t,iC} = summed_LFP_single{t,iC};
        mean_signal=mean(summed_LFP_cleaned{t,iC},2);
        std_signal=std(mean_signal);
        idxToRemove=find(abs(mean_signal)>std_signal*2);
        summed_LFP_cleaned{t,iC}(idxToRemove,:)=[];
        lfpSTA{t,iC} = mean(summed_LFP_cleaned{t,iC});
    end
end

%% Channels To Select
channelsToPlot = [58];
ind = 660:710;
ampArray = {};
brainArea = 'M1';

% 1 = individual, 2 = tiled, 3 = heatmap, 4 = multi-trial 5 = peak to peak
saveFormat = 5;

if brainArea == "S1"
% Array 1072-1
    
    mapping = [33 45 57 27 15 3 24 12
               35 47 59 25 13 1 22 10
               37 49 61 23 11 32 20 8
               39 51 63 21 9 30 18 6
               41 53 31 19 7 28 16 4
               43 55 29 17 5 26 14 2];
    channels = s1channels;
else
% Array 1072-14
   
    mapping = [66 78 90 69 81 93 42 54
               68 80 92 71 83 95 44 56
               70 82 94 73 85 34 46 58
               72 84 96 75 87 36 48 60
               74 86 65 77 89 38 50 62
               76 88 67 79 91 40 52 64];
    channels = m1channels;
end
    

for t=1:length(TrialsToLoad)
    time = [-pre:1/Fs*1000:post];
    
     for aa = 1:length(channelsToPlot) 
         iC = (channelsToPlot(aa));

         singleSTA = lfpSTA{t,iC};
         
        if saveFormat == 1
            figure;
            plot(singleSTA)
            ylabel('Amplitude (uV)')
            xlabel('Time (mSec)')

            if sum(iC==s1channels) == 1
                area = 's1';
            else
                area = 'm1';
            end

            title(string(area) +' PAD '+string(TrialsToLoad(t))+ ' Electrode '+' ('+string(channelsToLoad(iC))+')'+' LFP)')
       
        end
     end


  if saveFormat == 2
        figure;
        count=1;
        axes = [];
        for row = 1:6
            for col = 1:8
                ax = subplot(6,8,count);
                ch = mapping(row, col);
                single_ch = lfpSTA{t,ch};
                hold on
                plot(time,single_ch)
                
                
                count = count+1;
                if brainArea == 'S1'
                    title(string(brainArea) +' PAD '+string(TrialsToLoad(t))+ ' Elec '+' ('+string(channelsToLoad(ch))+')'+' STA')
        
                else
                    title(string(brainArea) +' PAD '+string(TrialsToLoad(t))+ ' Elec '+' ('+string(channelsToLoad(ch))+')'+' STA')
                
                end
                axes(end+1) = ax;
            end
        end
        linkaxes(axes,'xy')
  end

     if saveFormat==3
        amplitudes = plotHeatMap(brainArea,TrialsToLoad(t),ind,lfpSTA(t,:),mapping);
        figure;
        boxplot(amplitudes)
        histogram(amplitudes,10)
        file = 'PAD_'+string(TrialsToLoad(t))+'_'+string(brainArea) +'_EP_Amps.mat';
        save(file,'amplitudes')
%         
     end

    if saveFormat == 5
        
        %calculate peak to peak amplitude
        for iC = 1:length(channels)
            single_ch = summed_LFP_cleaned{t,iC};
            for aa = 1:size(single_ch,1)
                win = single_ch(aa,670:1300); % IC antidromic window
                amp = peak2peak(win);
                ampArray{t,iC}(aa)=amp;
            end
        end
   end

end
if saveFormat ==5
file = 'PAD_'+string(TrialsToLoad(1))+'_'+string(brainArea) +'_EP_Amps.mat';
    save(file,'ampArray')    
end

if saveFormat == 4
    figure;
    axes = [];
   [x,y] = optSubplotLayout(length(channelsToPlot)*2);
   for t = 1:length(TrialsToLoad)
       if t == 1
           c= [.7 .7 .7];
           cAvg = [0 0 0];
           shift = 0;
       else
           c= [.7 .7 .7];
           cAvg = [0 0 0];
           shift = 1;
       end
           for aa = 1:length(channelsToPlot)
               iC = channelsToPlot(aa);
               ax = subplot(x,y,aa+shift);
               hold on
               title('PAD '+string(TrialsToLoad(1))+' '+string(brainArea) + ' Elec '+' ('+string(iC)+')'+' STA')
               
               if t==2
                   plot(time,summed_LFP_cleaned{t,iC}(200:250,:),'color',c);
               else
                   plot(time,summed_LFP_cleaned{t,iC}(20:60,:),'color',c);
               end
               plot(time,lfpSTA{t,iC},'color',cAvg,'LineWidth',1.5);
               
           
               axes(end+1) = ax;
               xlim([-5 40])
               ylim([-3000 6000])
           end
%           legend
   end
   linkaxes(axes,'xy')

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

function vector=plotHeatMap(brainArea,trial,ind, all_LFP_mat,mapping)
vector = [];
for row = 1:6
   for col = 1:8
        
        ch = mapping(row, col);
        single_ch = all_LFP_mat{ch};
        win = single_ch(624:1100);
        p2p = peak2peak(win);

        amp = max(single_ch(ind));
        latMat(row,col) = amp;

        arrayMatrix(row,col) = p2p;
        vector(end+1) = p2p;
       
   end
end

custom = zeros(48,3);
left_color = [1 1 1]; % white
right_color = [.87 0 0]; % color according to the attached image
redcol = interp1([0, 1], [left_color; right_color], linspace(0, 1, 45));

greycol = linspace(.5,1,3);

custom(1:3,1) = greycol;
custom(1:3,2) = greycol;
custom(1:3,3) = greycol;
custom(4:48,:) = redcol;
% custom(:,:) = redcol;


figure;
imagesc(arrayMatrix,[500 7000])
colormap(custom)
colorbar
colormapeditor
axis square;

title(string(brainArea) +' PAD '+string(trial))

colormap(custom)
colorbar
colormapeditor
axis square;

title(string(brainArea) +' PAD '+string(trial))

figure;
heatmap(arrayMatrix)

end

function filtered_LFP = notchFilter_LFP(LFP)
fs_chan = 30000;

filterbands_Line=[58,62];
[z, p, k] = butter(2, filterbands_Line/(fs_chan/2), 'stop');
[sos_line,g_line] = zp2sos(z, p, k, 'down', 'two');

filterbands=[10,5000];
[z,p,k] = butter(2, filterbands/(fs_chan/2), 'bandpass');
[sos,g] = zp2sos(z, p, k, 'down', 'two');


filtered_LFP=filtfilt(sos,g,double(LFP'))';
filtered_LFP=filtfilt(sos_line,g_line,double(filtered_LFP'))';
end   

function [xNum,yNum] = optSubplotLayout(numChn)
% Determine the most equal layout for the subplot figures

% Determine the number of subplots
xFig = 16;
yFig = 9;

%b = factor(numChn);
cond = 1:numChn./2;%b(1):b(end);
area = nan(length(cond),1);

for n = 1:length(cond)
    temp = [cond(n) ceil(numChn/cond(n))];
    area(n,1) = abs((xFig/max(temp))-(yFig/min(temp)));
    area(n,2) = [temp(1)*temp(2)-numChn]./numChn;
    area(n,3) = (max(temp)-min(temp))./numChn;
    pairDat(n,:) = [max(temp) min(temp)];
    area(n,4) = .05*area(n,1) + .475*(area(n,2) + area(n,3));
end

idx1 = find(area(:,4) == min(area(:,4)));

if length(idx1)>1
    [val,idx2] = min(area(idx1,1));
    idx = idx1(idx2);
else
    idx = idx1;
end

xNum = pairDat(idx,1);
yNum = pairDat(idx,2);
end


function [dat_ns5] = loadExperimentData(FILENAME,PATHNAME)
% Automatically loads in all the raw data from the ripple data for the
% acute experiment.

% Sampling rates
sr_ns5 = 30000; % Analog data/broadband
sr_ns2 = 1000; % LFP
sr_nf3 = 2000;

aD = 'ns5'; % Identify the analog data channels

% Load the analog data
dat_ns5 = openNSxCervical(fullfile(PATHNAME,[FILENAME '.' aD]));


end