%load('C:\Users\emg27\Dropbox\PostDoc\ElviraMarco\VOP_Paper\figure4\panelC_EMGData_20230801.mat')
load('D:\Dropbox\PostDoc\ElviraMarco\VOP_Paper\figure4\panelC_EMGData_20230801.mat')
% Calculate the area under the curve plot for each muscle and condition
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
figC2(1) = figure;
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
figC2(2) = figure;
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
figC2(3) = figure;
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
