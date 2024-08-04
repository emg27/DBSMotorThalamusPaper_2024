clear  
close all

%% Open FILE
pathName = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-JC M68-22 Johnny Cage\20221005 - Experiment\DATA\mat\';
TrialsToLoad = [3 12 13 14 15 17 18 19 20];
% TrialsToLoad = [3];
% channelsToLoad = [1:128];
% TrialsToLoad = [3 12 8 10];
channelsToLoad = [183 189 164 179	187	185	168	175	191	181	170	165	162	177	180	161	166	173	188	169	172	163	174	167	176	171	182	186	184	190	192	178	195	201	205	222	199	197	209	218	203	193	219	216	207	224	223	206	211	220	215	198	221	214	217	212	213	210	200	204	196	202	208	194];

analog = 2; % analog ch 2 (IC STIM) or 3 (VOP STIM)
pre = 5;
post = 50;

% Array 1072-1
    
    S1mapping = [97 103 117 122 109 128 9 11
               98 104 118 106 125 112 13 15
               99 113 119  123 110 1 2 17
               100 114 120 107 126 3 4 16
               101 115 105 108 111 5 6 10
               102 116 121 124 127 7 8 12];
    S1channels = reshape(S1mapping,1,[]);

    % Array 1072-14
   
    M1mapping = [19 26 49 41 44 48 65 69
               21 25 33 40 45 57 66 77
               18 30 34 52 43 59 73 70
               14 28 50 38 55 60 74 78
               24 27 36 39 56 61 67 71
               23 32 35 42 58 62 75 79
               22 29 51 53 46 63 68 72
               20 31 37 54 47 64 76 80];
     M1channels = reshape(M1mapping,1,[]);


%% grabs timestamps from NS5 and LFP from 
StimThresh=16000;
Fs = 30000;

for t=1:length(TrialsToLoad)
    ns5_filename = sprintf('NS5_PAD%04d',TrialsToLoad(t));
    load([pathName ns5_filename '.mat']);
    
   stimAll{t}=findAnalogStim(NS5,StimThresh);
   StimTimes{t}=stimAll{t}(stimAll{t}(:,1)==analog,2);

    for iC = 1:length(channelsToLoad) 
        ch = channelsToLoad(iC);
        
        rawLFP{t,iC} = double(NS5.Data(ch,:));
        filteredLFP{t,iC} = notchFilter_LFP(rawLFP{t,iC});
    end   
end

    %% Build STA
time = -pre:1/30:post;
for t=1:length(TrialsToLoad)
    for iC = 1:length(channelsToLoad)
        singleLFP = filteredLFP{t,iC};
        for aa = 1:length(StimTimes{t})
           stim = StimTimes{t}(aa);
           lowerbound = stim-(pre*Fs/1000);
           upperbound = stim+(post*Fs/1000);

           if upperbound < length(singleLFP) || lowerbound < 0
               window = singleLFP(lowerbound:upperbound);
               %baseline correct
               window = window - mean(window(1:(stim-lowerbound-100))); 
%                window(594:667) = 0; % 1ms blanking of stim artifact, IC antidromic

               summed_LFP_single{t,iC}(aa,:) = (window);

               allP2P{t,iC}(aa,:) = peak2peak(window(270:600));
           end
        end 
    end   

    for iC = 1:length(channelsToLoad)
    %Clean Variables
        summed_LFP_cleaned{t,iC} = summed_LFP_single{t,iC};
        mean_signal=mean(summed_LFP_cleaned{t,iC},2);
        std_signal=std(mean_signal);
        idxToRemove=find(abs(mean_signal)>std_signal);
%         summed_LFP_cleaned{t,iC}(idxToRemove,:)=[];
        lfpSTA{t,iC} = mean(summed_LFP_cleaned{t,iC});
    end
end

        

%% BOXPLOTS
close all
% delay = {'IC Alone'}
delay = {'IC Alone','2ms','5ms','10ms','50ms','100ms','150ms','300ms','500ms'};
% delay = {'IC Alone','2ms','50 Hz','100Hz'};
channelsToPlot = [183 189 164 179 187	185	168	175	191	181	170	165	162	177	180	161	166	173	188	169	172	163	174	167	176	171	182	186	184	190	192	178 195	201	205	222	199	197	209	218	203	193	219	216	207	224	223	206	211	220	215	198	221	214	217	212	213	210	200	204	196	202	208	194];

for ii = 1:length(channelsToPlot)
    figure;
    axes = [];
    p = [];
    p2pMat = [];
    ind = find(channelsToLoad==channelsToPlot(ii));

    for t=1:length(TrialsToLoad)
        singleP2P = allP2P{t,ind};
        p2pCell{t} = singleP2P;      
        p2pMat = [p2pMat; p2pCell{t}];
        p = [p; ones(size(p2pCell{t}))*t];
    
    end
    
     subplot(1,1,1)
    boxplot(p2pMat,p,'Labels',delay,'symbol','+r')
    title('Ch '+string(channelsToPlot(ii))+' Spine P2P')
    
    trialMat = [1 2
                ];
    
    alpha = .001;
    offset = 0;
    for comb = 1:size(trialMat,1)
        trial1 = trialMat(comb,1);
        trial2 = trialMat(comb,2);
    
        [CI,sig] = bootstrapCompMeans(allP2P{trial1,ind},allP2P{trial2,ind},10000,alpha/size(trialMat,1));
    
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
    saveas(gca,'Johnny_conditioningDelays_spineEp_ch'+string(channelsToPlot(ii)))
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
    
%     Plot histogram to confirm
%     figure;
%     hist(diffSampMeans, 100)
    
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
