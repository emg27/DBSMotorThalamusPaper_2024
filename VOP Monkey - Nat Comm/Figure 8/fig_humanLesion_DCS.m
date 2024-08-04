% Figure 8 Plot: Demonstrate that the results are consistent in the
% impaired patient
clear
close all

dataset = '20220405';
saveDir = uigetdir;
setNum = {'#11','#12'}; %{'#11','#12'};
useAvg = 1;
useMA = 1;
numSmth = 21;
spacePR = 20;
colS1 = [81 163 191]/252;
colM1 = [36 62 140]/252;
idxChnECOG = [2 4];
offsetECOG = 32;
PATHNAME = 'D:\Figures\VOP_NatCom_2024_Data\XLTEK';%
%% PANEL Phase Revisal plots - supplemental

contacts = {'[Fz - 1]','[Fz - 2]','[Fz - 3]','[Fz - 4]','[FZ - 5]',...
    '[FZ - 6]','[FZ - P3]'};

% Collect the reference contact
refDat = XLTEK.export(dataset,PATHNAME,'testType','Phase Reversal',...
    'contact',contacts{end},'setNum',setNum);

% Collect the ECOG contacts
Fpr = figure; hold on
for n = 1:6
    data{n} = XLTEK.export(dataset,PATHNAME,'testType','Phase Reversal',...
        'contact',contacts{n},'setNum',setNum);
    temp = refDat-data{n};
    % Determine if we are averaging across plots
    if useAvg
        temp = mean(temp,2);
    end
    % Determine if we are using smoothing
    if useMA % Ie use box car
        temp = movmean(temp,numSmth);
    end
    temp = temp - temp(1,:); % Align zero to the first point
    plot(linspace(0,100,600),temp+(6-n)*spacePR,'k')
end
line([80 80],[0 10],'Color','r')
text(81,6.5,'10 \muV')
% Correct the labels
xlabel('Time (ms)')
ylabel('Contact Number')
gx = gca;
set(gx, 'YTick',[0:5]*spacePR,'YTickLabel',{'1','2','3','4','5','6'});
title(sprintf('Subject %s, Phase Reversal',dataset))
Fpr.Name = sprintf('%s_phaseReversal',dataset);

%% PANEL B - Create the ECOG alignment plots
chnECOG = offsetECOG + idxChnECOG;
chnID = {'ECOG 1','ECOG 2','ECOG 3','ECOG 4','ECOG 5','ECOG 6'};
PLTTITLE = {'2 Hz Stim'};
filename = 'datafile0007';
fs = 30000;
stimThresh = 4000;
fbands = [10 500];%[10 2000];
numEx = 30;

% Timing Information
preStim = 10;%5;
postStim = 20;%30;
tBlnk = 3; 
dt = 1000/fs;
tStep = -preStim:dt:postStim;
tBOff = 0.7;
binW = [.25]; % Rank sum one tail side.
bsVal = tBlnk:.5:postStim;

% Load the data and find the stmulation times
[dat_ns5] = loadExperimentData(filename,'P:\projects\human\VOP STIM\DATA\Intra_Op_DBS\20220405');%fullfile(PATHNAME,dataset)); %'
[stim] = findAnalogStim(sum(dat_ns5.Data(offsetECOG+[1:6],:)),stimThresh);
if isempty(stim)
    stim = [1:fs/2:size(dat_ns5.Data,2)]';
    stim = [0*stim stim];
end

%% Proces and plot the data
checkFilt.dat = ECOG.calcSTA(dat_ns5.Data(chnECOG,:),chnID(idxChnECOG),stim,...
    'fbands',fbands,'plotData',0,'pre_t',preStim,'post_t',postStim);

% Run rank sum for the data
[pRS hRS] = deal(nan(length(bsVal),1));

for n = 1:length(bsVal)
    idx = tStep>=(bsVal(n)-binW) & tStep<(bsVal(n)+binW);
    data1 = double(checkFilt.dat{1}(:,idx==1)) - mean(checkFilt.dat{1}(:,tStep<-tBOff),2);
    data2 = double(checkFilt.dat{2}(:,idx==1)) - mean(checkFilt.dat{2}(:,tStep<-tBOff),2);
    [pRS(n),hRS(n)] = ranksum(data1(:),data2(:));%,'tail','left'); %06/29 edit
end

% Plot the data
idx = randperm(size(stim,1),numEx);
s1Trl = double(checkFilt.dat{1})-mean(checkFilt.dat{1}(:,tStep<-tBOff),2);
m1Trl = double(checkFilt.dat{2})-mean(checkFilt.dat{2}(:,tStep<-tBOff),2);
Fecog = figure;
Fecog.Name = sprintf('%s_10Hz_ecogComp',dataset);
hold on
pltS1Trl = plot(tStep,s1Trl(idx,:),'Color',[colS1 0.3]);
pltM1Trl = plot(tStep,m1Trl(idx,:),'Color',[colM1 0.3]);

pltS1Avg = plot(tStep,mean(s1Trl),'Color',colS1,'LineWidth',2);
pltM1Avg = plot(tStep,mean(m1Trl),'Color',colM1,'LineWidth',2);

% Label the figure
title('10Hz Comparison')
xAxLim = [-preStim postStim];
axis([xAxLim -250 125])

% Add the significance to the plot
figure(Fecog)
sigplt = plot(bsVal(hRS==1),100,'kd');
xline(-tBOff,'k')
xline(tBlnk,'k')

legend([pltM1Avg pltS1Avg sigplt(1)],'M1 ECOG','S1 ECOG','p<0.05','Location','best')

% Plot the ECOGs separately
Fecog(2) = figure;
Fecog(2).Name = sprintf('%s_10Hz_ecogComp_Separate',dataset);

% S1
subplot(2,1,1)
hold on
plot(tStep,s1Trl(idx,:),'Color',[colS1 0.3]);
plot(tStep,mean(s1Trl),'Color',colS1,'LineWidth',2);

% M1
subplot(2,1,2)
hold on
plot(tStep,m1Trl(idx,:),'Color',[colM1 0.3]);
plot(tStep,mean(m1Trl),'Color',colM1,'LineWidth',2);
plot(bsVal(hRS==1),100,'kd');
xline(-tBOff,'k')
xline(tBlnk,'k')

saveFigurePDF(Fecog,saveDir)

%% PANEL B: Peak to peak box plot of the motor thalamus stimulation only
s1Trunk = s1Trl(:,tStep>tBlnk);
m1Trunk = m1Trl(:,tStep>tBlnk);

% Calculate
maxS1 = max(s1Trunk,[],2);
maxM1 = max(m1Trunk,[],2);
minS1 = min(s1Trunk,[],2);
minM1 = min(m1Trunk,[],2);
p2pS1 = maxS1 - minS1;
p2pM1 = maxM1 - minM1;

% Remove outliers
[~,idxS1Out] = rmoutliers(p2pS1);
[~,idxM1Out] = rmoutliers(p2pM1);
p2pS1(idxS1Out==1) = nan;
p2pM1(idxM1Out==1) = nan;

[p_p2p,h_p2p] = ranksum(p2pS1,p2pM1);%,'tail','left'); %06/29 edit

% Plot the boxplot
F = figure;
F.Name = 'ECOG_p2p_boxplot';
boxplot([p2pS1 p2pM1])
title(sprintf('p = %d, h = %f0.0 ',p_p2p,h_p2p))
saveFigurePDF(F,saveDir)
%% PANEL C: Plot the DCS figures
contactAll = {'APB','FLEX'};%,'EXT'}; % Best are APB, FLEX, EXT
allCond = {'13mA','12mA','11mA'};
PATHFILE = 'D:\Dropbox\PostDoc\ElviraMarco\IntraOperativeForms\XLTEK';%'C:\Users\emg27\Dropbox\PostDoc\ElviraMarco\IntraOperativeForms\XLTEK\';%'C:\Users\emg27\Dropbox\PostDoc\ElviraMarco\IntraOperativeForms\XLTEK\';%
allCond = {'13mA','12mA','11mA'};
time = linspace(0,100,600);
dt = mean(diff(time));
xLimPlt = [20 75];
aocAll = repmat(struct('Contact',[],'Stim',[],'dM1',[],'d50',[],...
    'd100',[],'p',[],'h',[]),length(contactAll),length(allCond));
%saveDir = uigetdir;

for n = 1:length(contactAll)
    contact = contactAll(n)
    for k = 3%1:length(allCond)
        cond = allCond{k};
        switch cond
            case '13mA'
                % Contact 3: VOP DBS Stim Off 13mA (383-390)
                setM1 = 383:390;
                % Contact 3, VOP 4V 50Hz 100us, 13mA 422-445
                set50 = 422:445;
                % COntact 3, VOP: 4V 100Hz 100us, 13mA 409-421
                set100 = 409:421;
            case '12mA'
                % Contact 3: VOP DBS Stim Off 12mA (327-346) - no rmv cannules
                setM1 = 327:346;
                % Contact 3, VOP 4V 50Hz 100us, 12mA 446-458
                set50 = 446:458;
                % COntact 3, VOP: 4V 100Hz 100us, 12mA 510-529
                set100 = 510:529;
            case '11mA'
                % Contact 3: VOP DBS Stim Off 11mA (364-382)
                setM1 = 364:382;
                % Contact 3, VOP 4V 50Hz 100us, 11mA: 478-486
                set50 = 478:486;
                % COntact 3, VOP: 4V 100Hz 100us, 11mA 530-549
                set100 = 530:549;
        end
        lM1= length(setM1);
        l50 = length(set50);
        l100 = length(set100);

        [aocM1, aoc50, aoc100] = deal(nan(1,max([lM1 l50 l100])));

        % 10 mA M1 only
        dataM1 = XLTEK.export(dataset,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',setM1'));

        % VOP DBS On (3mA 50 Hz) 10 mA M1
        data50 = XLTEK.export(dataset,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',set50'));

        % VOP DBS On (3mA 100 Hz) 10 mA M1
        data100 = XLTEK.export(dataset,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',set100'));

        Fdcs(k,1) = figure('Position',[100 100 800 850]);
        axeD(1) = subplot(3,2,1); hold on
        plot(time,dataM1)
        plot(time,mean(dataM1,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n M1 stim only %s',contact{1},cond))

        axeR(1) = subplot(3,2,2); hold on
        plot(time,abs(dataM1))
        plot(time,mean(abs(dataM1),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n M1 stim only %s',contact{1},cond))

        axeD(2) = subplot(3,2,3); hold on
        plot(time,data50)
        plot(time,mean(data50,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n VOP 50Hz + M1 %s',contact{1},cond))

        axeR(2) = subplot(3,2,4); hold on
        plot(time,abs(data50))
        plot(time,mean(abs(data50),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n VOP 50Hz + M1 %s',contact{1},cond))

        axeD(3) = subplot(3,2,5); hold on
        plot(time,data100)
        plot(time,mean(data100,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n VOP 100Hz + M1 %s',contact{1},cond))
        xlabel('Index'), ylabel('Trace')

        axeR(3) = subplot(3,2,6); hold on
        plot(time,abs(data100))
        plot(time,mean(abs(data100),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n VOP 100Hz + M1 %s',contact{1},cond))

        matchAxis(axeD)
        setLimits.xLim = xLimPlt;
        matchAxis(axeD,'setLimits',setLimits,'setY',0,'setZ',0,'setR',0,'setTheta',0)
        matchAxis(axeR)
        matchAxis(axeR,'setLimits',setLimits,'setY',0,'setZ',0,'setR',0,'setTheta',0)
        Fdcs(k,1).Name = sprintf('%s_dcs_traces_%s_m1_%s',dataset,contact{1},cond);
        Fdcs(k,1).Name = sprintf('%s_dcs_traces_%s_m1_%s',dataset,contact{1},cond);

        % Subselect the indices to include in the area under the curve
        subset = 125:400;%1:600;%150:300;

        aocM1(1:lM1) = trapz(time(subset),abs(dataM1(subset,:)));
        aoc50(1:l50) = trapz(time(subset),abs(data50(subset,:)));
        aoc100(1:l100) = trapz(time(subset),abs(data100(subset,:)));
        
        % Store the data
        aocAll(n,k).Contact = contact{:};
        aocAll(n,k).Stim = cond;
        aocAll(n,k).dM1 = aocM1;
        aocAll(n,k).d50 = aoc50;
        aocAll(n,k).d100 = aoc100;
        
        % Run rank sum statistics
        [p(1),h(1)] = ranksum(aocAll(n,k).dM1,aocAll(n,k).d50);
        [p(2),h(2)] = ranksum(aocAll(n,k).dM1,aocAll(n,k).d100);
        [p(3),h(3)] = ranksum(aocAll(n,k).d50,aocAll(n,k).d100);
        aocAll(n,k).p = p;
        aocAll(n,k).h = h;
        
        Fdcs(k,2) = figure('Position',[900 100 800 850]);
        boxplot([aocM1; aoc50; aoc100]')
        title(sprintf('Area under the curve: %s \n Stim Settings: M1 %s',contact{1},cond))
        xlabel('Conditions'), ylabel('Area under the curve')
        set(gca,'XTickLabels',{'M1 Only','VOP 50 Hz','VOP 100 Hz'})
        Fdcs(k,2).Name = sprintf('%s_dcs_auc_%s_m1_%s',dataset,contact{1},cond);

        % Plot trace response through time
        Fdcs(k,3) = figure;%('Position',[100 100 800 850]);
        axeT(1) = subplot(3,1,1);
        d = dataM1(subset,:);
        plot(0:dt:dt*[length(d(:))-1],d(:))
        %plot(abs(d(:)))
        xlabel('Time (ms)'), ylabel('Trace')
        title(sprintf('Rect. DCS through time %s \n M1 stim only %s',contact{1},cond))

        axeT(2) = subplot(3,1,2);
        d = data50(subset,:);
        plot(0:dt:dt*[length(d(:))-1],d(:))
        %plot(abs(d(:)))
        xlabel('Time (ms)'), ylabel('Trace')
        title(sprintf('Rect. DCS through time %s \n VOP 50Hz + M1 %s',contact{1},cond))

        axeT(3) = subplot(3,1,3);
        d = data100(subset,:);
        plot(0:dt:dt*[length(d(:))-1],d(:))
        %plot(abs(d(:)))
        xlabel('Time (ms)'), ylabel('Trace')
        title(sprintf('Rect. DCS through time %s \n VOP 100Hz + M1 %s',contact{1},cond))
        
        matchAxis(axeT)
        Fdcs(k,3).Name = sprintf('%s_dcs_time_%s_m1_%s',dataset,contact{1},cond);
    end
end