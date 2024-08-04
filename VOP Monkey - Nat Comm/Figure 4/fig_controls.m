%% Figure 4 Plots: 1) radial EMG vs Radial + VOP emg, p2p values auc values,
% 3) EMG: VOP + IC + Post IC (lesion-pre post)

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

%% PANEL C: EMG comparisons of Radial Nerve and VOP

% Load the data - OPAL
[datOpal,FOpal] = calEvokedResponse(pad_o,cond_o,'OPAL',PATHNAME_o,...
    'stimIdx',3,'sortedch',sortedchOpal,'EMGoffset',0,'idxoffset',128,...
    'bankoffset',128,'calcSpine',0);
datOpal(6).condition = cond_o(end); 

% Combine pad 22 and pad 65
[datOpal] = combinePAD(datOpal,5,6);
[datOpal] = combinePAD(datOpal,[2 3],5);

% Load the data - Hosu
[datHosu,FHosu] = calEvokedResponse(pad_h,cond_h,'HOSU',PATHNAME_h,...
    'stimIdx',2,'sortedch',sortedchHosu,'EMGoffset',96,'idxoffset',0,...
    'bankoffset',0,'calcSpine',0);

% Load the data - Johnny
[datJohnny,FJohnny] = calEvokedResponse(pad_j,cond_j,'JOHNNY',PATHNAME_j,...
    'stimIdx',2,'sortedch',sortedchJohnny,'EMGoffset',128,'idxoffset',0,...
    'bankoffset',256-160,'calcSpine',0);
% Combine pad 5 and pad 168
[datJohnny] = combinePAD(datJohnny,[1 2],7);
% Combine pad 25 and pad 26
[datJohnny] = combinePAD(datJohnny,[5 6],8);

close([FOpal,FHosu,FJohnny])
%% Panel C1: Create the example traces
emgChn_O = [1:10];
condIdx_O = [1 5 4 6]; % 10Hz, Radial, Radial+VOP, 50Hz
stimIdx = [2 3 3 2];
colMat = [VOPpair; radial; VOPpair; VOPpair];

% Opal Plots
[F_emgVOP(1)] = plotRadialComp(datOpal,condIdx_O,emgChn_O,stimIdx,colMat,preStim,postStim)
F_emgVOP(1).Name = 'radialResponse_Opal';
%ylim([-50 50])

% Panel C1: Create the Hosu example trace
emgChn_H = [1 2 4 5 7:10];
condIdx_H = [1 3 4 5]; % 10Hz, 50Hz, Radial, Radial+VOP
stimIdx = [3 3 2 2];
colMat = [VOPpair; VOPpair; radial; VOPpair];

% Opal Plots
[F_emgVOP(3)] = plotRadialComp(datHosu,condIdx_H,emgChn_H,stimIdx,colMat,preStim,postStim)
F_emgVOP(3).Name = 'radialResponse_Hosu';
%ylim([-50 50])

% Panel C1: Create the Hosu example trace
emgChn_J = [1 2 4 5 6 7:10];
condIdx_J = [2 3 4 8]; % 10Hz, 50Hz, Radial, Radial+VOP
stimIdx = [3 3 2 2];
colMat = [VOPpair; VOPpair; radial; VOPpair];

% Opal Plots
[F_emgVOP(2)] = plotRadialComp(datJohnny,condIdx_J,emgChn_J,stimIdx,colMat,preStim,postStim)
F_emgVOP(2).Name = 'radialResponse_Johnny';
%ylim([-50 50])

%% Calculate the force peak to peak for each muscle and condition
% Iterate through each muscle of each animal

% Opal
emgChn_O = [4 3];
%emgChn_O = [3:5 7];%1:10;%
%condIdx_O =[1 5 4]; % 10Hz, Radial, Radial+VOP
[p2pVal_O,grp_O] = deal(cell(length(emgChn_O),length(condIdx_O)));
for n = 1:length(emgChn_O)
    for k = 1:length(condIdx_O)
        p2pVal_O{n,k} = calPkToPk(datOpal(condIdx_O(k)).emg(emgChn_O(n)).indTrl)';
        grp_O{n,k} = k*ones(size(p2pVal_O{n,k}));
    end
    
    % Run the bootstrap stats between the conditions,
    [ciBoot_O{n,1},rnBoot_O(n,1)] = bootstrapCompMeans(p2pVal_O{n,1},p2pVal_O{n,2},10000,.001,1,'l');
    [ciBoot_O{n,2},rnBoot_O(n,2)] = bootstrapCompMeans(p2pVal_O{n,1},p2pVal_O{n,3},10000,.001,1,'l');
    [ciBoot_O{n,3},rnBoot_O(n,3)] = bootstrapCompMeans(p2pVal_O{n,2},p2pVal_O{n,3},10000,.001,1,'l');
    
    % Run the Wilcoxon Rand Sum stats
    [p_O(n,1),h_O(n,1)] = ranksum(p2pVal_O{n,1},p2pVal_O{n,2},'tail','left');
    [p_O(n,2),h_O(n,2)] = ranksum(p2pVal_O{n,1},p2pVal_O{n,3},'tail','left');
    [p_O(n,3),h_O(n,3)] = ranksum(p2pVal_O{n,2},p2pVal_O{n,3},'tail','left');
end
% Plot the data
figC2(1) = figure;
figC2(1).Name = 'Opal_radialNerve_p2p';
for n = 1:length(emgChn_O)
    subplot(1,length(emgChn_O),n)
    boxplot([p2pVal_O{n,:}],[grp_O{n,:}])
    title(['EMG ' num2str(emgChn_O(n))])
end
plotTitle('Opal')


% Hosu
emgChn_H = [2 1];
%emgChn_H = [1 2 4 5 7];%1:10;%
%condIdx_H =[1 3 4 5]; % 10Hz, 50Hz, Radial, Radial+VOP
[p2pVal_H,grp_H] = deal(cell(length(emgChn_H),length(condIdx_H)));
for n = 1:length(emgChn_H)
    for k = 1:length(condIdx_H)
        p2pVal_H{n,k} = calPkToPk(datHosu(condIdx_H(k)).emg(emgChn_H(n)).indTrl)';
        grp_H{n,k} = k*ones(size(p2pVal_H{n,k}));
    end
    
    % Run the bootstrap stats between the conditions,
    [ciBoot_H{n,1},rnBoot_H(n,1)] = bootstrapCompMeans(p2pVal_H{n,1},p2pVal_H{n,2},10000,.001,1,'l');
    [ciBoot_H{n,2},rnBoot_H(n,2)] = bootstrapCompMeans(p2pVal_H{n,2},p2pVal_H{n,3},10000,.001,1,'l');
    [ciBoot_H{n,3},rnBoot_H(n,3)] = bootstrapCompMeans(p2pVal_H{n,2},p2pVal_H{n,4},10000,.001,1,'l');
    [ciBoot_H{n,4},rnBoot_H(n,4)] = bootstrapCompMeans(p2pVal_H{n,3},p2pVal_H{n,4},10000,.001,1,'l');
    
    % Run the Wilcoxon Rand Sum stats
    [p_H(n,1),h_H(n,1)] = ranksum(p2pVal_H{n,1},p2pVal_H{n,2},'tail','left');
    [p_H(n,2),h_H(n,2)] = ranksum(p2pVal_H{n,2},p2pVal_H{n,3},'tail','left');
    [p_H(n,3),h_H(n,3)] = ranksum(p2pVal_H{n,2},p2pVal_H{n,4},'tail','left');
    [p_H(n,4),h_H(n,4)] = ranksum(p2pVal_H{n,3},p2pVal_H{n,4},'tail','left');
end
% Plot the data
figC2(2) = figure;
figC2(2).Name = 'Hosu_radialNerve_p2p';
for n = 1:length(emgChn_H)
    subplot(1,length(emgChn_H),n)
    boxplot([p2pVal_H{n,:}],[grp_H{n,:}])
    title(['EMG ' num2str(emgChn_H(n))])
end
plotTitle('Hosu')

% Johnny
emgChn_J = [1 2];
%emgChn_J = [1 2 4 5 6 7];%1:10;%
%condIdx_J =[2 3 4 8]; % 10Hz, 50Hz, Radial, Radial+VOP
[p2pVal_J,grp_J] = deal(cell(length(emgChn_J),length(condIdx_J)));
for n = 1:length(emgChn_J)
    for k = 1:length(condIdx_J)
        p2pVal_J{n,k} = calPkToPk(datJohnny(condIdx_J(k)).emg(emgChn_J(n)).indTrl)';
        grp_J{n,k} = k*ones(size(p2pVal_J{n,k}));
    end
    
    % Run the bootstrap stats between the conditions,
    [ciBoot_J{n,1},rnBoot_J(n,1)] = bootstrapCompMeans(p2pVal_J{n,1},p2pVal_J{n,2},10000,.001,1,'l');
    [ciBoot_J{n,2},rnBoot_J(n,2)] = bootstrapCompMeans(p2pVal_J{n,2},p2pVal_J{n,3},10000,.001,1,'l');
    [ciBoot_J{n,3},rnBoot_J(n,3)] = bootstrapCompMeans(p2pVal_J{n,2},p2pVal_J{n,4},10000,.001,1,'l');
    [ciBoot_J{n,4},rnBoot_J(n,4)] = bootstrapCompMeans(p2pVal_J{n,3},p2pVal_J{n,4},10000,.001,1,'l');
    
    % Run the Wilcoxon Rand Sum stats
    [p_J(n,1),h_J(n,1)] = ranksum(p2pVal_J{n,1},p2pVal_J{n,2},'tail','left');
    [p_J(n,2),h_J(n,2)] = ranksum(p2pVal_J{n,2},p2pVal_J{n,3},'tail','left');
    [p_J(n,3),h_J(n,3)] = ranksum(p2pVal_J{n,2},p2pVal_J{n,4},'tail','left');
    [p_J(n,4),h_J(n,4)] = ranksum(p2pVal_J{n,3},p2pVal_J{n,4},'tail','left');
end
% Plot the data
figC2(3) = figure;
figC2(3).Name = 'Johnny_radialNerve_p2p';
for n = 1:length(emgChn_J)
    subplot(1,length(emgChn_J),n)
    boxplot([p2pVal_J{n,:}],[grp_J{n,:}])
    title(['EMG ' num2str(emgChn_J(n))])
end
plotTitle('Johnny')

%saveFigurePDF([figC2(:); F_emgVOP(:)])

%% Plot the VOP responses only

figE1(1) = figure;
figE1(1).Name = 'Opal_VOPOnly_compare_p2p';
for n = 1:length(emgChn_O)
    tmpGrp = [grp_O{n,:}];
    tmpP2P = [p2pVal_O{n,:}];
       idx = ismember(tmpGrp,[1 4]);
    subplot(1,length(emgChn_O),n)
    boxplot(tmpP2P(idx==1),tmpGrp(idx==1))
    title(['EMG ' num2str(emgChn_O(n))])
    ylim([0 45])
end
plotTitle('Opal')

figE1(2) = figure;
figE1(2).Name = 'Johnny_VOPOnly_compare_p2p';
for n = 1:length(emgChn_J)
    tmpGrp = [grp_J{n,:}];
    tmpP2P = [p2pVal_J{n,:}];
       idx = ismember(tmpGrp,[1 2]);
    subplot(1,length(emgChn_J),n)
    boxplot(tmpP2P(idx==1),tmpGrp(idx==1))
    title(['EMG ' num2str(emgChn_J(n))])
    ylim([0 45])
end
plotTitle('Johnny')

figE1(3) = figure;
figE1(3).Name = 'Hosu_VOPOnly_compare_p2p';
for n = 1:length(emgChn_H)
    tmpGrp = [grp_H{n,:}];
    tmpP2P = [p2pVal_H{n,:}];
       idx = ismember(tmpGrp,[1 2]);
    subplot(1,length(emgChn_H),n)
    boxplot(tmpP2P(idx==1),tmpGrp(idx==1))
    title(['EMG ' num2str(emgChn_H(n))])
    ylim([0 45])
end
plotTitle('Hosu')

%% Plot the VOP and radial responses only

figE2(1) = figure;
figE2(1).Name = 'Opal_RadialVOPOnly_compare_p2p';
for n = 1:length(emgChn_O)
    tmpGrp = [grp_O{n,:}];
    tmpP2P = [p2pVal_O{n,:}];
       idx = ismember(tmpGrp,[2 3]);
    subplot(1,length(emgChn_O),n)
    boxplot(tmpP2P(idx==1),tmpGrp(idx==1))
    title(['EMG ' num2str(emgChn_O(n))])
end
plotTitle('Opal')

figE2(2) = figure;
figE2(2).Name = 'Johnny_RadialVOPOnly_compare_p2p';
for n = 1:length(emgChn_J)
    tmpGrp = [grp_J{n,:}];
    tmpP2P = [p2pVal_J{n,:}];
       idx = ismember(tmpGrp,[3 4]);
    subplot(1,length(emgChn_J),n)
    boxplot(tmpP2P(idx==1),tmpGrp(idx==1))
    title(['EMG ' num2str(emgChn_J(n))])
end

figE2(3) = figure;
figE2(3).Name = 'Hosu_RadialVOPOnly_compare_p2p';
for n = 1:length(emgChn_H)
    tmpGrp = [grp_H{n,:}];
    tmpP2P = [p2pVal_H{n,:}];
       idx = ismember(tmpGrp,[3 4]);
    subplot(1,length(emgChn_H),n)
    boxplot(tmpP2P(idx==1),tmpGrp(idx==1))
    title(['EMG ' num2str(emgChn_H(n))])
end
plotTitle('Hosu')
%% Extended figure: Rerun the bootstrap with different alpha values
alpha = [0.05 0.01 0.001]

exFigDat = repmat(struct('alpha',[],'ci_o',[],'h_o',[],'ci_h',[],'h_h',[],...
    'ci_j',[],'h_j',[]),length(alpha),1);
for k = 1:length(alpha)
    a = alpha(k);
    
    % Opal
    emgChn_O = [4 3];%[3:5 7];%1:10;%
    condIdx_O =[1 5 4 6]; % 10Hz, Radial, Radial+VOP
    for n = 1:length(emgChn_O)
        % Run the bootstrap stats between the conditions,
        [ciBoot_O{n,1},rnBoot_O(n,1)] = bootstrapCompMeans(p2pVal_O{n,1},p2pVal_O{n,2},10000,a,2,'l');
        [ciBoot_O{n,2},rnBoot_O(n,2)] = bootstrapCompMeans(p2pVal_O{n,1},p2pVal_O{n,3},10000,a,2,'l');
        [ciBoot_O{n,3},rnBoot_O(n,3)] = bootstrapCompMeans(p2pVal_O{n,2},p2pVal_O{n,3},10000,a,2,'l');
        [ciBoot_O{n,4},rnBoot_O(n,4)] = bootstrapCompMeans(p2pVal_O{n,1},p2pVal_O{n,4},10000,a,2,'l');
    end
    
    % Hosu
    emgChn_H = [2 1];
    %emgChn_H = [1 2 4 5 7]; %1:10;%
    condIdx_H =[1 3 4 5]; % 10Hz, 50Hz, Radial, Radial+VOP
    for n = 1:length(emgChn_H)
        % Run the bootstrap stats between the conditions,
        [ciBoot_H{n,1},rnBoot_H(n,1)] = bootstrapCompMeans(p2pVal_H{n,1},p2pVal_H{n,2},10000,a,2,'l');
        [ciBoot_H{n,2},rnBoot_H(n,2)] = bootstrapCompMeans(p2pVal_H{n,2},p2pVal_H{n,3},10000,a,2,'l');
        [ciBoot_H{n,3},rnBoot_H(n,3)] = bootstrapCompMeans(p2pVal_H{n,2},p2pVal_H{n,4},10000,a,2,'l');
        [ciBoot_H{n,4},rnBoot_H(n,4)] = bootstrapCompMeans(p2pVal_H{n,3},p2pVal_H{n,4},10000,a,2,'l');
    end
    
    % Johnny
    emgChn_J = [1 2];
    %emgChn_J = [1 2 4 5 6 7]; %1:10;%
    condIdx_J =[2 3 4 8]; % 10Hz, 50Hz, Radial, Radial+VOP
    for n = 1:length(emgChn_J)
        % Run the bootstrap stats between the conditions,
        [ciBoot_J{n,1},rnBoot_J(n,1)] = bootstrapCompMeans(p2pVal_J{n,1},p2pVal_J{n,2},10000,a,2,'l');
        [ciBoot_J{n,2},rnBoot_J(n,2)] = bootstrapCompMeans(p2pVal_J{n,2},p2pVal_J{n,3},10000,a,2,'l');
        [ciBoot_J{n,3},rnBoot_J(n,3)] = bootstrapCompMeans(p2pVal_J{n,2},p2pVal_J{n,4},10000,a,2,'l');
        [ciBoot_J{n,4},rnBoot_J(n,4)] = bootstrapCompMeans(p2pVal_J{n,3},p2pVal_J{n,4},10000,a,2,'l');
    end
    
    % Fill the structure
    exFigDat(k).alpha = a;
    exFigDat(k).ci_o = ciBoot_O;
    exFigDat(k).h_o = rnBoot_O;
    exFigDat(k).ci_h = ciBoot_H;
    exFigDat(k).h_h = rnBoot_H;
    exFigDat(k).ci_j = ciBoot_J;
    exFigDat(k).h_j = rnBoot_J;
end

%% Panel D: Calculate the area under the curve plot for each muscle and condition
% Interate through each muscle of each animal
alpha = [0.05 0.01 0.001];

preStim = 5; %30; % ms
postStim = 15;%100; % ms
time = -preStim:1/30:postStim;
minT = 1;
maxT = 10;
idx = find(time>minT & time<=maxT);

% Opal
emgChn_O = [4 3];
%emgChn_O = [3:5 7]; %1:10;%
condIdx_O = [1 5 4]; % 10Hz, Radial. Radial+VOP
[aucVal_O,grp_O] = deal(cell(length(emgChn_O),length(condIdx_O)));
for n = 1:length(emgChn_O)
    for k = 1:length(condIdx_O)
        temp = abs(datOpal(condIdx_O(k)).emg(emgChn_O(n)).indTrl);
        aucVal_O{n,k} = trapz(time(idx),temp(:,idx)');
        grp_O{n,k} = k*ones(size(aucVal_O{n,k}));
    end

    % Run the bootstrap stats between the conditions,
    for a = 1:length(alpha)
        [ciBoot_Oa{n,1,a},rnBoot_Oa(n,1,a)] = bootstrapCompMeans(aucVal_O{n,1},aucVal_O{n,2},10000,alpha(a),1);
        [ciBoot_Oa{n,2,a},rnBoot_Oa(n,2,a)] = bootstrapCompMeans(aucVal_O{n,1},aucVal_O{n,3},10000,alpha(a),1);
        [ciBoot_Oa{n,3,a},rnBoot_Oa(n,3,a)] = bootstrapCompMeans(aucVal_O{n,2},aucVal_O{n,3},10000,alpha(a),1);
    end
end
% Plot the data
fig2C2(1) = figure;
fig2C2(1).Name = 'Opal_radialNerve_auc';
for n = 1:length(emgChn_O)
    subplot(1,length(emgChn_O),n)
    boxplot([aucVal_O{n,:}],[grp_O{n,:}])
    title(['EMG ' num2str(emgChn_O(n))])
end
plotTitle('Opal')

% Hosu
emgChn_H = [2 1];
%emgChn_H = [1 2 4 5 7]; %1:10;%
condIdx_H = [1 3 4 5]; % 10Hz, 50Hz, Radial. Radial+VOP
[aucVal_H,grp_H] = deal(cell(length(emgChn_H),length(condIdx_H)));
for n = 1:length(emgChn_H)
    for k = 1:length(condIdx_H)
        temp = abs(datHosu(condIdx_H(k)).emg(emgChn_H(n)).indTrl);
        aucVal_H{n,k} = trapz(time(idx),temp(:,idx)');
        grp_H{n,k} = k*ones(size(aucVal_H{n,k}));
    end

    % Run the bootstrap stats between the conditions,
    for a = 1:length(alpha)
        [ciBoot_Ha{n,1,a},rnBoot_Ha(n,1,a)] = bootstrapCompMeans(aucVal_H{n,1},aucVal_H{n,2},10000,alpha(a),1);
        [ciBoot_Ha{n,2,a},rnBoot_Ha(n,2,a)] = bootstrapCompMeans(aucVal_H{n,2},aucVal_H{n,3},10000,alpha(a),1);
        [ciBoot_Ha{n,3,a},rnBoot_Ha(n,3,a)] = bootstrapCompMeans(aucVal_H{n,2},aucVal_H{n,4},10000,alpha(a),1);
        [ciBoot_Ha{n,4,a},rnBoot_Ha(n,4,a)] = bootstrapCompMeans(aucVal_H{n,3},aucVal_H{n,3},10000,alpha(a),1);
    end
end
% Plot the data
fig2C2(2) = figure;
fig2C2(2).Name = 'Hosu_radialNerve_auc';
for n = 1:length(emgChn_H)
    subplot(1,length(emgChn_H),n)
    boxplot([aucVal_H{n,:}],[grp_H{n,:}])
    title(['EMG ' num2str(emgChn_H(n))])
end
plotTitle('Hosu')

% Johnny
emgChn_J = [1 2];
%emgChn_J = [1 2 4 5 6 7]; %1:10;%
condIdx_J = [2 3 4 8]; % 10Hz, 50Hz, Radial. Radial+VOP
[aucVal_J,grp_J] = deal(cell(length(emgChn_J),length(condIdx_J)));
for n = 1:length(emgChn_J)
    for k = 1:length(condIdx_J)
        temp = abs(datJohnny(condIdx_J(k)).emg(emgChn_J(n)).indTrl);
        aucVal_J{n,k} = trapz(time(idx),temp(:,idx)');
        grp_J{n,k} = k*ones(size(aucVal_J{n,k}));
    end

    % Run the bootstrap stats between the conditions,
    for a = 1:length(alpha)
        [ciBoot_Ja{n,1,a},rnBoot_Ja(n,1,a)] = bootstrapCompMeans(aucVal_J{n,1},aucVal_J{n,2},10000,alpha(a),1);
        [ciBoot_Ja{n,2,a},rnBoot_Ja(n,2,a)] = bootstrapCompMeans(aucVal_J{n,2},aucVal_J{n,3},10000,alpha(a),1);
        [ciBoot_Ja{n,3,a},rnBoot_Ja(n,3,a)] = bootstrapCompMeans(aucVal_J{n,2},aucVal_J{n,4},10000,alpha(a),1);
        [ciBoot_Ja{n,4,a},rnBoot_Ja(n,4,a)] = bootstrapCompMeans(aucVal_J{n,3},aucVal_J{n,3},10000,alpha(a),1);
    end
end
% Plot the data
fig2C2(3) = figure;
fig2C2(3).Name = 'Johnny_radialNerve_auc';
for n = 1:length(emgChn_J)
    subplot(1,length(emgChn_J),n)
    boxplot([aucVal_J{n,:}],[grp_J{n,:}])
    title(['EMG ' num2str(emgChn_J(n))])
end
plotTitle('Johnny')

%saveFigurePDF([figc2(:); F_emgVOP(:)])

%% Extended figure 4: Calculate the area under the curve plot for each 
% muscle at 10 and 50Hz conditions
% Interate through each muscle of each animal
alpha = [0.05 0.01 0.001];

preStim = 5; %30; % ms
postStim = 15;%100; % ms
time = -preStim:1/30:postStim;
minT = 1;
maxT = 10;
idx = find(time>minT & time<=maxT);

% Opal
emgChn_O = [1:10];
condIdx_O = [1 6]; % 10Hz, 50Hz
[aucVal_O,grp_O] = deal(cell(length(emgChn_O),length(condIdx_O)));
for n = 1:length(emgChn_O)
    for k = 1:length(condIdx_O)
        temp = abs(datOpal(condIdx_O(k)).emg(emgChn_O(n)).indTrl);
        aucVal_O{n,k} = trapz(time(idx),temp(:,idx)');
        grp_O{n,k} = k*ones(size(aucVal_O{n,k}));
    end

    % Run the bootstrap stats between the conditions,
    for a = 1:length(alpha)
        [ciBoot_Oa{n,1,a},rnBoot_Oa(n,1,a)] = bootstrapCompMeans(aucVal_O{n,1},aucVal_O{n,2},10000,alpha(a),1);
    end
end
% Plot the data
figEx4(1) = figure;
figEx4(1).Name = 'Opal_VOPOnly_compare_auc';
for n = 1:length(emgChn_O)
    subplot(1,length(emgChn_O),n)
    boxplot([aucVal_O{n,:}],[grp_O{n,:}])%, 'DataLim',[-inf 250])
    title(['EMG ' num2str(emgChn_O(n))])
    ylim([0 200])
end
plotTitle('Opal')

% Hosu
emgChn_H = [1 2 4 5 7:8];
condIdx_H = [1 3]; % 10Hz, 50Hz
[aucVal_H,grp_H] = deal(cell(length(emgChn_H),length(condIdx_H)));
for n = 1:length(emgChn_H)
    for k = 1:length(condIdx_H)
        temp = abs(datHosu(condIdx_H(k)).emg(emgChn_H(n)).indTrl);
        aucVal_H{n,k} = trapz(time(idx),temp(:,idx)');
        grp_H{n,k} = k*ones(size(aucVal_H{n,k}));
    end

    % Run the bootstrap stats between the conditions,
    for a = 1:length(alpha)
        [ciBoot_Ha{n,1,a},rnBoot_Ha(n,1,a)] = bootstrapCompMeans(aucVal_H{n,1},aucVal_H{n,2},10000,alpha(a),1);
    end
end
% Plot the data
figEx4(2) = figure;
figEx4(2).Name = 'Hosu_VOPOnly_compare_auc';
for n = 1:length(emgChn_H)
    subplot(1,length(emgChn_H),n)
    boxplot([aucVal_H{n,:}],[grp_H{n,:}])%, 'DataLim',[-inf 250])
    title(['EMG ' num2str(emgChn_H(n))])
    ylim([0 200])
end
plotTitle('Hosu')

% Johnny
emgChn_J = [1 2 4 5 6 7 9];
condIdx_J = [2 3]; % 10Hz, 50Hz
[aucVal_J,grp_J] = deal(cell(length(emgChn_J),length(condIdx_J)));
for n = 1:length(emgChn_J)
    for k = 1:length(condIdx_J)
        temp = abs(datJohnny(condIdx_J(k)).emg(emgChn_J(n)).indTrl);
        aucVal_J{n,k} = trapz(time(idx),temp(:,idx)');
        grp_J{n,k} = k*ones(size(aucVal_J{n,k}));
    end

    % Run the bootstrap stats between the conditions,
    for a = 1:length(alpha)
        [ciBoot_Ja{n,1,a},rnBoot_Ja(n,1,a)] = bootstrapCompMeans(aucVal_J{n,1},aucVal_J{n,2},10000,alpha(a),1);
    end
end
% Plot the data
figEx4(3) = figure;
figEx4(3).Name = 'Johnny_VOPOnly_compare_auc';
for n = 1:length(emgChn_J)
    subplot(1,length(emgChn_J),n)
    boxplot([aucVal_J{n,:}],[grp_J{n,:}])%, 'DataLim',[-inf 250])
    title(['EMG ' num2str(emgChn_J(n))])
    ylim([0 200])
end
plotTitle('Johnny')

%saveFigurePDF([figEx4(:)])

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

function [p2pVal,maxVal,minVal] = calPkToPk(a,varargin)
preStim = 5; % ms
postStim = 15; % ms
t_blnk = 0.5; % ms

assignopts(who, varargin);

t = -preStim:1/30:postStim;

idx = t>=t_blnk;
a = a(:,idx==1);
maxVal = max(a,[],2);
minVal = min(a,[],2);
p2pVal = maxVal - minVal;

% Remove Outliers
[~,idx_rm] = rmoutliers(p2pVal);
p2pVal(idx_rm == 1) = nan;
maxVal(idx_rm == 1) = nan;
minVal(idx_rm == 1) = nan;
end

%% Create EMG example plots
function [F_emgVOP] = plotRadialComp(dat,condIdx,emgChn,stimIdx,colMat,preStim,postStim)
col = length(condIdx);
row = length(emgChn);

F_emgVOP = figure('Position',[14 98 1008 897]);
for k = 1:col
    pos = find(dat(condIdx(k)).stimAll(:,1) == stimIdx(k));
    if length(pos)>30
        idx = randperm(length(pos),30);
    else
        idx = 1:length(pos)
    end
    
    % EMG examples: FOr Opal use emgs 4 and 5
    figure(F_emgVOP)
    for n = 1:length(emgChn)
        ax(n,k) = subplot(row,col,k+(n-1)*col), hold on
        plot(-preStim:1/30:postStim,dat(condIdx(k)).emg(emgChn(n)).indTrl(idx,:),...
            'Color',[colMat(k,:) 0.18],'LineWidth',.25)
        plot(-preStim:1/30:postStim,nanmean(dat(condIdx(k)).emg(emgChn(n)).indTrl),....
            'Color',colMat(k,:),'LineWidth',1)
        if n == 1
            title(dat(condIdx(k)).condition)
        end
        xlim([-preStim postStim])
        xlabel('Time (ms)'), ylabel(['EMG ' num2str(emgChn(n))])
        %if n == row
        linkaxes(ax(n,:))
                
        %end
    end
end
end