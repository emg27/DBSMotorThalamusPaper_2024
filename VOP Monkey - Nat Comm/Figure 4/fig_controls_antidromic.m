%% Figure 4 Plots: 1) Antidromic response in the spinal cord no emg signal
% VOP alone, 2) radial EMG vs Radial + VOP emg, p2p values, 3) EMG: VOP +
% IC + Post IC (lesion-pre post)

% Default settings
sr_ns5 = 30000;
preStim = 5;%30; % ms
postStim = 15;%100; % ms
stimThresh = 16000; % Stimulation detection

% Colors
VOPpair = [36 62 140]/252;
radial = [57 181 74]/252;

% PADs to analyze in figure 4
PATHNAME_o = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-OP M63-18 Opal\20211208 - Experiments\DATA';
pad_o = [21 22 65 47 11];
cond_o = ["VOP 10Hz, 21","RADIAL, 22","RADIAL, 65","VOP + RADIAL, 47", "VOP 50 Hz, 11"];
sortedchOpal = [18 16 20 14 17 15 19 13 31 1 29 3 32 2 27 5 30 4 25 7 28 6 23 9 26 8 21 11 24 10 22 12];

PATHNAME_h = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-HS M219-21 Hosu\20220427 - Experiment\DATA';
pad_h = [21 28 23 26 27];
cond_h = ["VOP 10Hz  PAD 21","VOP 10Hz  PAD 28",...
    "VOP 50Hz  PAD 23", "RADIAL  PAD 26","VOP + RADIAL  PAD 27"];
% probe 1 (right) without jumper
sortedch1 = 128+[10	4 29 14	6 8	25 18 2	12 23 28 31	16 13 32 27	20 5 24	21 30 19 26	17 22 11 7	9 3	1 15];
% probe 2 (left) with jumper
sortedch2 = 160+[30 24 20 3 26 28 16 7 22 32 6 9 18 1 2 19 14 5 10 27 4 11 8 13 12 15 25 21 29 23 17 31];
sortedchHosu = [sortedch1 sortedch2];

PATHNAME_j = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-JC M68-22 Johnny Cage\20221005 - Experiment\DATA';
pad_j = [5 168 163 24 25 26];
cond_j = ["VOP 10Hz  PAD 5","VOP 10Hz  PAD 168","VOP 50Hz  PAD 163",...
    "RADIAL  PAD 24","VOP + RADIAL  PAD 25","VOP + RADIAL  PAD 26"];
% Spinal cord antidromic
sortedch1 = 160+[23,29,4,19,27,25,8,15,31,21,10,5,2,17,20,1,6,13,28,9,12,3,14,7,16,11,22,26,24,30,32,18];
% probe 2 (left) with jumper
sortedch2 = 160+32+[3,9,13,30,7,5,17,26,11,1,27,24,15,32,31,14,19,28,23,6,29,22,25,20,21,18,8,12,4,10,16,2];
sortedchJohnny = [sortedch1 sortedch2];

%% Load the data for panels A and B
[datOpal,FOpal] = calEvokedResponse(21,cond_o(1),'OPAL',PATHNAME_o,...'stimIdx',2,
    'sortedch',sortedchOpal,'EMGoffset',0,'idxoffset',128,'bankoffset',128);
[datHosu,FHosu] = calEvokedResponse(21,cond_h(1:2),'HOSU',PATHNAME_h,...'stimIdx',3,
    'sortedch',sortedchHosu,'EMGoffset',96,'idxoffset',0,'bankoffset',0);
[datJohnny,FJohnny] = calEvokedResponse([5 168],cond_j(1:2),'JOHNNY',...'stimIdx',3,
    PATHNAME_j,'sortedch',sortedchJohnny,'EMGoffset',128,'idxoffset',0,'bankoffset',256-160);
close([FOpal,FHosu,FJohnny])

%% PANEL A: antidromic plot examples from Opal
chns = [1 16 19 31];
pos = find(datOpal(1).stimSTA(:,1) == 2);
idx = randperm(length(pos),30);

% Antidromic Response: 10 Hz VOP alone
F_spineVOP = figure('Position',[165 50 318 913]);
F_spineVOP.Name = 'spine_antidromicResponse_Opal';
figure(F_spineVOP)
for n = 1:length(chns)
    subplot(length(chns),1,n),hold on
    plot(-preStim:1/30:postStim,datOpal(1).spine{chns(n)}(pos(idx),:),...
        'Color',[[36 62 140]/252 0.18],'LineWidth',.25)
    plot(-preStim:1/30:postStim,nanmean(datOpal(1).spine{chns(n)}(pos,:)),....
        'Color',[36 62 140]/252,'LineWidth',1)
    if n == 1
        title(sprintf('%s \n Channel %0.0f (D to V)',datOpal(1).condition,chns(n)))
    else
        title(sprintf('Channel %0.0f (D to V)',chns(n)))
    end
    xlim([-2 10])
    xlabel('Time (ms)'),ylabel('Spine')
    axis square
    set(gca,'XTick',-2:10)
end
linkaxes

%% PANEL B: Quantification of antidromic response, all animals
% CalcPeaktoPeak data
[amp_o,width_o,tP_o,tN_o] = calPkToPkSpine(datOpal(1),2,0,7);
[ampB_o,widthB_o,tPB_o,tNB_o] = calPkToPkSpine(datOpal(1),2,0,7);

[amp_h,width_h,tP_h,tN_h] = calPkToPkSpine(datHosu(1),3,1,7);

[amp_j1,width_j1,tP_j1,tN_j1] = calPkToPkSpine(datJohnny(1),3,1,7);
[amp_j2,width_j2,tP_j2,tN_j2] = calPkToPkSpine(datJohnny(2),3,1,7);
amp_j = [amp_j1 amp_j2];
width_j = [width_j1 width_j2];
%%
tP_j = [tP_j2];% tP_j2];
tN_j = [tN_j2];% tN_j2];

% Combine the peak to peak
p2pVal_O = [nanmean(amp_o,2) nanmean(amp_j(1:32,:),2) nanmean(amp_h(1:32,:),2)];
p2pVal_O(9,1) = nan;

% Create the colormap
colMat = ones(length(unique(p2pVal_O(:))),3);
colMat(:,2) = linspace(1,0,length(unique(p2pVal_O(:))));
colMat(:,3) = linspace(1,0,length(unique(p2pVal_O(:))));

%% Plot the heatmap separated by animal
titleVal = {'Opal','Johnny','Hosu'};
F(2) = figure;
F(2).Name = 'allAnimals_spinalResponse_Peak2PeakAmp_separate';
for n = 1:3
    subplot(1,3,n)
    colormap(colMat)
    heatmap(p2pVal_O(:,n))
    colormap(colMat)
    
    %imagesc(p2pVal_O(:,n))
    % h = colorbar;
    %h.Label.String = 'Normalized Amplitude'
    ylabel('Channel Number')
    title(titleVal{n})
end

figure; heatmap((p2pVal_O-min(p2pVal_O))./(max(p2pVal_O)-min(p2pVal_O)))
colormap(colMat)
F(6) = gcf
F(6).Name = 'allAnimals_spinalResponse_Peak2PeakAmp_separate_v4';

%% Plot the time latency
F(3) = figure;
F(3).Name = 'allAnimals_antidromic_latency'
boxplot([nanmean(tP_o,2) nanmean(tP_j(1:32,:),2) nanmean(tP_h(1:32,:),2)])
ylim([0 5])
ylabel('Antidromic Onset (ms)')
set(gca,'XTick',1:3,'XTickLabel',{'Opal','Johnny','Hosu'},'XTickLabelRotation',30)

%saveFigurePDF(F)

%% Plot the time latency
F(3) = figure;
F(3).Name = 'allAnimals_antidromic_latency'
myboxplot({nanmean(tP_o,1)'; nanmean(tP_j(1:32,:),1)'; nanmean(tP_h(1:32,:),1)'},'box')
ylim([0 5])
ylabel('Antidromic Onset (ms)')
set(gca,'XTick',1:3,'XTickLabel',{'Opal','Johnny','Hosu'},'XTickLabelRotation',30)
%% Functions

function [dat,F] = calEvokedResponse(PAD,conditions,subject,PATHNAME,varargin)
stimIdx = 2;
sortedch = [];
sr_ns5 = 30000;
preStim = 5;%30; % ms
postStim = 15;%100; % ms
stimThresh = 16000; % Stimulation detection
EMGoffset = 0;
idxoffset = 160;
bankoffset = 256;
optLoadArg = [];
baseTrue = 1;
fBands = [10 5000];
calcSpine = 1;

assignopts(who, varargin);

if isempty(sortedch)
    sortedch = 1:64;
end

F = figure;
for n = 1:length(PAD)
    
    FILENAME = sprintf('PAD%04d',PAD(n));
    if ~isempty(optLoadArg)
        [dat_ns5,nev,~] = loadExperimentData(FILENAME,PATHNAME,optLoadArg);
    else
        [dat_ns5] = loadExperimentData(FILENAME,PATHNAME,'loadNev',0);
    end
    
    stim = findAnalogStim(dat_ns5,stimThresh);
    dat(n).subject = subject;
    dat(n).PAD = PAD(n);
    dat(n).condition = conditions(n);
    dat(n).stimAll = stim;
    
    % EMG Data
    if size(stim(stim(:,1)==stimIdx),1) == 0
        [dat(n).emg] = Hosu.spkTrigAvg(dat_ns5,stim(stim(:,1)~=stimIdx,:),...
            preStim,postStim,sr_ns5,EMGoffset + [1:20],F,...
            'emgChns',EMGoffset + [1:20],'baseTrue',baseTrue);
    else
        [dat(n).emg] = Hosu.spkTrigAvg(dat_ns5,stim(stim(:,1)==stimIdx,:),...
            preStim,postStim,sr_ns5,EMGoffset + [1:20],F,...
            'emgChns',EMGoffset + [1:20],'baseTrue',baseTrue);
    end
    
    
    if calcSpine
        try
            [dat(n).spine,~,~,~,dat(n).stimSTA] = ECOG.calcSTA(dat_ns5.Data(idxoffset+sortedch,:),...
                string(bankoffset+sortedch),stim,'optPlot',1,'fbands',fBands,...
                'post_t',postStim,'pre_t',preStim,'plotData',0);
        end
    end
end
end

function [dat] = combinePAD(dat,comIdx,finIdx)
dat(finIdx) = dat(comIdx(1));

dat(finIdx).subject = dat(comIdx).subject;
dat(finIdx).PAD = [dat(comIdx).PAD];
dat(finIdx).condition = sprintf('CombinedRadial PAD [%s]',num2str(dat(finIdx).PAD));

for k = 2:length(comIdx)
    dat(finIdx).stimAll = [dat(finIdx).stimAll; dat(comIdx(k)).stimAll];
    
    if isfield(dat,'spine')
        dat(finIdx).stimSTA = [dat(finIdx).stimSTA; dat(comIdx(k)).stimSTA];
        for n = 1:32
            dat(finIdx).spine{n} = [dat(finIdx).spine{n}; dat(comIdx(k)).spine{n}];
        end
    end
    
    for n = 1:10
        dat(finIdx).emg(n).indTrl = [dat(finIdx).emg(n).indTrl; dat(comIdx(k)).emg(n).indTrl];
        dat(finIdx).emg(n).avg = nanmean(dat(finIdx).emg(n).indTrl);
    end
end
end

%% Calculate the peak to peak values for Opal's VOP Only stim
function [ampTot,widthTot,tP,tN] = calPkToPkSpine(a,stimIdx,minT,maxT)
%condIdx =[2 9 5];% 7;
%stimIdx = [2 3 3];
preStim = 5;
postStim = 15;
%minT = 0;
%maxT = 7;

stim = find(a.stimSTA(:,1) == stimIdx(1));
x = -preStim:1/30:postStim;

[ampP,tP,widthP,ampN,tN,widthN] = deal(nan(size(a.spine,2),size(stim,1)));
for chn = 1:size(a.spine,2)
    for trl = 1:size(stim,1)
        % Calculate the positive peak information
        [ampPd,tPd,widthPd] = findpeaks(double(a.spine{chn}(stim(trl),x>minT & x<=maxT)),x(x>minT & x<=maxT),'SortStr','Descend');
        ampP(chn,trl) = ampPd(1);
        tP(chn,trl) = tPd(1);
        widthP(chn,trl) = widthPd(1);
        
        % Calculate the negative peak information
        [ampNd,tNd,widthNd] = findpeaks(-double(a.spine{chn}(stim(trl),x>minT & x<=maxT)),x(x>minT & x<=maxT),'SortStr','Descend');
        ampN(chn,trl) = ampNd(1);
        tN(chn,trl) = tNd(1);
        widthN(chn,trl) = widthNd(1);
    end
end
ampTot = ampP+ampN;
widthTot = widthN + widthP;

% Remove Outliers
[~,idx_rm] = rmoutliers(ampTot);
ampTot(idx_rm == 1) = nan;
widthTot(idx_rm == 1) = nan;
tP(idx_rm == 1) = nan;
tN(idx_rm == 1) = nan;
end