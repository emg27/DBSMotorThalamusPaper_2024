clear, close all
%% Load the finalized dat file and make the histogram figure
load('D:\Figures\VOP_NatCom_2024_Data\Figure 5\SOURCE_FIG5_frequencyResponseData_allAnimals.mat')

%% Collect the data to create the histograms
pltCond = {'Attenuation','Suppression','Alternation','Initial','No modulation'};
conditions = {{'attentuation','attenuation'},'suppression','alternation',{'initial','intial'},...
    {'little response','no response'}};

uniVOPFreq = unique({dat.vopFreq});

%% Calculate the number of response conditions organized by VOP frequencies,
% All animals combined, Opal, and Hosu
[fsVOPCond] = srtCondData(dat,uniVOPFreq,conditions,1);
idxO = ismember([dat.subject],'Opal');
[fsVOPCond_O] = srtCondData(dat(idxO==1),uniVOPFreq,conditions,1);
idxH = ismember([dat.subject],'Hosu');
[fsVOPCond_H] = srtCondData(dat(idxH==1),uniVOPFreq,conditions,1);
idxJ = ismember([dat.subject],'Johnny');
[fsVOPCond_J] = srtCondData(dat(idxJ==1),uniVOPFreq,conditions,1);
%% Identify the muscle type indices for each animal
% OPAL
% EMG 1 - APL, EMG 2 - FDM, EMG 3 - EDC, EMG 4 - ECR, EMG 5 - FCR
% EMG 6 - BIC, EMG 7 - FDC, EMG 8 - Jaw, EMG 9 - Superior lip, EMG 10 - Cheek
o_APL = ismember([dat.EMG],1) & idxO; % HAND
o_FDM = ismember([dat.EMG],2) & idxO; % HAND
o_EDC = ismember([dat.EMG],3) & idxO; % HAND
o_ECR = ismember([dat.EMG],4) & idxO; % HAND
o_FCR = ismember([dat.EMG],5) & idxO; % HAND
o_BIC = ismember([dat.EMG],6) & idxO; % ARM
o_FDC = ismember([dat.EMG],7) & idxO; % HAND
o_Jaw = ismember([dat.EMG],8) & idxO; % FACE
o_supLip = ismember([dat.EMG],9) & idxO; % FACE
o_cheek = ismember([dat.EMG],10) & idxO; % FACE

% HOSU
% EMG 1 - EDC, EMG 2 - ECR, EMG Error - FCR, EMG 3 - FDC, EMG 4 - ADP
% EMG 5 - FDM, EMG 6 - BIC, EMG 7 - Superior Lip, EMG 8 - Cheek, EMG 9 - Jaw
h_EDC = ismember([dat.EMG],1) & idxH; % HAND
h_ECR = ismember([dat.EMG],2) & idxH; % HAND
h_FDC = ismember([dat.EMG],3) & idxH; % HAND
h_ADP = ismember([dat.EMG],4) & idxH; % HAND
h_FDM = ismember([dat.EMG],5) & idxH; % HAND
h_BIC = ismember([dat.EMG],6) & idxH; % ARM
h_supLip = ismember([dat.EMG],7) & idxH; % FACE
h_cheek = ismember([dat.EMG],8) & idxH; % FACE
h_Jaw = ismember([dat.EMG],9) & idxH; % FACE

% Johnny
% EMG 1 - ECR, EMG 2 - EDC, EMG 3 (err) - FDC, EMG 4 (err) - FCR, EMG 5 - ABP
% EMG 6 - FDM, EMG 7 - BIC, EMG 8 - Cheek, EMG 9 - Jaw, EMG 10 - Sup Lip.
j_ECR = ismember([dat.EMG],1) & idxJ; % HAND
j_EDC = ismember([dat.EMG],2) & idxJ; % HAND
j_FDC = ismember([dat.EMG],3) & idxJ; % HAND
j_FCR = ismember([dat.EMG],4) & idxJ; % HAND
j_ABP = ismember([dat.EMG],5) & idxJ; % HAND
j_FDM = ismember([dat.EMG],6) & idxJ; % HAND
j_BIC = ismember([dat.EMG],7) & idxJ; % ARM
j_cheek = ismember([dat.EMG],8) & idxJ; % FACE
j_Jaw = ismember([dat.EMG],9) & idxJ; % FACE
j_supLip = ismember([dat.EMG],10) & idxJ; % FACE


% Create the indices for different conditions
idxHand = o_APL | o_FDM | o_EDC | o_FDC | o_ECR | o_FCR |...
    h_EDC | h_FDC | h_ADP | h_FDM | h_ECR |...
    j_ECR | j_EDC | j_FDC | j_FCR | j_ABP | j_FDM;
idxArm =  o_BIC | h_BIC;
idxFace = o_Jaw | o_supLip | o_cheek |...
    h_supLip | h_cheek | h_Jaw |...
    j_supLip | j_cheek | j_Jaw;
idxFlex = o_FDM | o_FCR | o_FDC | h_FDC | h_FDM | j_FDC | j_FCR | j_FDM;
idxExt = o_EDC | o_ECR | h_EDC | h_ECR | j_ECR | j_EDC;

idxWrist = idxFlex | idxExt;
idxHand2 = idxHand & ~idxWrist;

% Collect the data
[fsVOPCond_Hand] = srtCondData(dat(idxHand==1),uniVOPFreq,conditions,1);
[fsVOPCond_Arm] = srtCondData(dat(idxArm==1),uniVOPFreq,conditions,1);
[fsVOPCond_Face] = srtCondData(dat(idxFace==1),uniVOPFreq,conditions,1);
[fsVOPCond_Flex] = srtCondData(dat(idxFlex==1),uniVOPFreq,conditions,1);
[fsVOPCond_Ext] = srtCondData(dat(idxExt==1),uniVOPFreq,conditions,1);
[fsVOPCond_Wrist] = srtCondData(dat(idxWrist==1),uniVOPFreq,conditions,1);
[fsVOPCond_Hand2] = srtCondData(dat(idxHand2==1),uniVOPFreq,conditions,1);

%% Collect the data to create the histograms - 20230623 - 3 condition example
pltCond = {'Attenuation + Suppression','Potentiation','No potentiation'};
conditions = {{'attentuation','attenuation','suppression'},...
    {'alternation','initial','intial','inital'},{'little response','no response','response'}};

uniVOPFreq = unique({dat.vopFreq});

% Calculate the number of response conditions organized by VOP frequencies,
% All animals combined, Opal, and Hosu
[fsVOPCond] = srtCondData(dat,uniVOPFreq,conditions,0);
X = categorical(pltCond);

[F(1)] = plotHistogram(fsVOPCond,X,'Assume potentiation for all uncertain cases (ie cases with multiple designations)- 3 conditions'); % All animals and conditions combined
F(1).Name = [F(1).Name '_allCombined_3conditions'];

%% Collect the data to create the histograms - 20230623 - 4 condition
pltCond = {'Attenuation','Suppression','Potentiation','No potentiation'};
conditions = {{'attentuation','attenuation'},'suppression',...
    {'alternation','initial','intial','inital'},{'little response','no response','response'}};

uniVOPFreq = unique({dat.vopFreq});

% Calculate the number of response conditions organized by VOP frequencies,
% All animals combined, Opal, and Hosu
[fsVOPCond] = srtCondData(dat,uniVOPFreq,conditions,0);
X = categorical(pltCond);

[F(1)] = plotHistogram(fsVOPCond,X,'Assume potentiation for all uncertain cases (ie cases with multiple designations) - 4 Conditions'); % All animals and conditions combined
F(1).Name = [F(1).Name '_allCombined_potentiationAssumed'];


%% Function Combine conditions
function [fsVOPCond_com] = combStimCond(fsVOPCond)
fsVOPCond_com(1,:) = sum(fsVOPCond(29:31,:)); % IC
fsVOPCond_com(2,:) = sum(fsVOPCond(1:2,:)); % 10 Hz
fsVOPCond_com(3,:) = sum(fsVOPCond(20:23,:)); % 50 Hz
fsVOPCond_com(4,:) = fsVOPCond(26,:); % 80 Hz
fsVOPCond_com(5,:) = sum(fsVOPCond([5 10:11],:)); % 100 Hz
fsVOPCond_com(6,:) = fsVOPCond(19,:); % 200 Hz
end

% Sort the data
function [fsVOPCond_com,fsVOPCond] = srtCondData(dat,uniVOPFreq,conditions,combineDat)

% Create a histogram organized by VOP condition
fsVOPCond = nan(size(uniVOPFreq,2),length(conditions));

for n = 1:size(uniVOPFreq,2)
    % Only select the relevant EMG
    mask = ismember({dat.vopFreq},uniVOPFreq(n));
    tempDat = dat(mask==1);
    for k = 1:length(conditions)
        % Collect all combined
        fsVOPCond(n,k) = sum(contains({tempDat.freqDes},conditions{k}));
    end
end

% Combine across stim conditions
if combineDat
    [fsVOPCond_com] = combStimCond(fsVOPCond);
else
    fsVOPCond_com(1,:) = sum(fsVOPCond(30:31,:)); % IC
    fsVOPCond_com(2,:) = fsVOPCond(1,:); % 10 Hz
    fsVOPCond_com(3,:) = fsVOPCond(20,:); % 50 Hz
    fsVOPCond_com(4,:) = fsVOPCond(26,:); % 80 Hz
    fsVOPCond_com(5,:) = fsVOPCond(5,:); % 100 Hz
    fsVOPCond_com(6,:) = fsVOPCond(19,:); % 200 Hzs
end
end

% Plot the histograms
function[F] = plotHistogram(fsVOPCond,X,TITLE)

F = figure;
F.Name = 'FrequencySuppressHistogram';
bar(X,100*fsVOPCond./sum(fsVOPCond,2))
legend('IC','10Hz VOP','50Hz VOP','80Hz VOP','100Hz VOP','200Hz VOP')
ylabel('Percentage of responses')
xlabel('Frequency Responses')
title(TITLE)
ylim([0 100])

end