%% Create the example traces of the data for Hosu and Opal
colors = [32 36 84; 35 64 138; 68 88 138; 99 111 138]./256;
postT = 10; %10
exPad_H = 24;
%filename = sprintf('arrayData_Hosu_n5top15_PAD%04d.mat',exPad);

if exist('C:\Users\emg27\Dropbox\')
    PATHNAME = 'C:\Users\emg27\Dropbox\PostDoc\ElviraMarco\VOP_Paper\figure5\arrayData\12-21\';
    PATHNAME_H = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-HS M219-21 Hosu\20220427 - Experiment\DATA\mat';
    PATHNAME_O = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-Opal M63-18 Opal\20211208 - Experiments\DATA\mat';
else
    PATHNAME = 'D:\Dropbox\PostDoc\ElviraMarco\VOP_Paper\figure5\arrayData\12-21\';
    PATHNAME_H = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-HS M219-21 Hosu\20220427 - Experiment\DATA\mat';
    PATHNAME_O = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-OP M63-18 Opal\20211208 - Experiments\DATA\mat';
end
F = figure; hold on
for monk = 1:2
    % Choose the animal information
    if monk == 1
        subject = 'Hosu';
        mapping_array1_m1_Hosu = [66 78 90 69 81 93	42 54;...
            68	80	92	71	83	95	44	56;...
            70	82	94	73	85	34	46	58;...
            72	84	96	75	87	36	48	60;...
            74	86	65	77	89	38	50	62;...
            76	88	67	79	91	40	52	64];
        chnM1 = reshape(mapping_array1_m1_Hosu,1,[]);
        
        chnEx = 42;
        pos = 1;%find(chnM1==chnEx);
        PATHNAME = PATHNAME_H;
        exPad = [21 23 24 25];
        cond = {'10Hz','50Hz','80Hz','100Hz'};
    else
        subject = 'Opal';
        mapping_array1_m1_Opal = [98 110 122 101 113 125 74	86;...
            100	112	124	103	115	127	76	88;...
            102	114	126	105	117	66	78	90;...
            104	116	128	107	119	68	80	92;...
            106	118	97	109	121	70	82	94;...
            108	120	99	111	123	72	84	96];
        
        chnM1 = reshape(mapping_array1_m1_Opal,1,[]);
        trlCnt = 80;
        chnEx = 78;%80;%99;%120;
        pos = 1;%find(chnM1==chnEx);
        PATHNAME = PATHNAME_O;
        exPad = [9 11 12 13];
        cond = {'10Hz','50Hz','80Hz','100Hz'};
    end
    
    showTrl = 20;
    time = -2:1/30:postT;
    
    
    for n = 1:length(exPad)
        % Load the data
        filename = sprintf('NS5_PAD%04d.mat',exPad(n));
        load(fullfile(PATHNAME,filename))
        
        % Run calc STA on just the one channel
        % Settings
        stimThresh = 8000;
        fBands = [20 5000];
        preStim = 2;%30; % ms
        postStim = 10;%100; % ms
        stimOn = findAnalogStim(NS5,stimThresh);
        stimOff = findAnalogStim(NS5,stimThresh,-1);
        [dat.brain,~,~,~,dat.stimSTA] = ECOG.calcSTA(NS5.Data(chnEx,:),...
            string(chnEx),stimOn,'fbands',fBands,'post_t',postT,'pre_t',2,...
            'plotData',0,'stimOff',stimOff,'flip',-1,'baseCor',1,'bw_n',2,'useSmooth',1);
       
        % Plot the first 80 trials
        trlCnt = size(dat.stimSTA,1);
        idxTrl = randperm(trlCnt,showTrl);
        
        % Plot the random trials
        figure(F)
        ax(monk,n) = subplot(2,4,n+4*(monk-1));
        hold on
        plot(time,dat.brain{pos}(idxTrl,:),'Color',[colors(n,:) .25])
        avg = mean(dat.brain{pos});
        plot(time,avg,'Color',colors(n,:),'LineWidth',2)
        yline(min(avg(time>=1.5)),'Color',colors(n,:),'LineWidth',2,'LineStyle','--')
        yline(max(avg(time>=1.5)),'Color',colors(n,:),'LineWidth',2,'LineStyle','--')
        xlim([-2 postT])
        title(sprintf('%s %s',subject,filename))

        % Calculate the auc of the data (For the single neuron)
        idx = find(time>=2);
        auc{monk,n} = trapz(time(idx),double(dat.brain{1}(:,idx))');
        aucIdx{monk,n} = n*ones(size(auc{monk,n}));

        % Save the trace data
        cep(monk,n).subject = subject;
        cep(monk,n).condition = cond{n};
        cep(monk,n).ch = chnEx;
        cep(monk,n).data = dat.brain{1,1};
        cep(monk,n).time = time;
    end
    matchAxis(ax(monk,:))

    figure
    b = [auc{monk,:}];
    b2 = [aucIdx{monk,:}];
    boxplot(b,b2)
    title(sprintf('%s %s AUC',subject,filename))
end

%% Create the peak to peak plots for each animal
for monk = 1:2
    % Choose the animal information
    if monk == 1
        subject = 'Hosu';
        mapping_array1_m1_Hosu = [66 78 90 69 81 93	42 54;...
            68	80	92	71	83	95	44	56;...
            70	82	94	73	85	34	46	58;...
            72	84	96	75	87	36	48	60;...
            74	86	65	77	89	38	50	62;...
            76	88	67	79	91	40	52	64];
        chnM1 = reshape(mapping_array1_m1_Hosu,1,[]);
        load('D:\Figures\VOP_NatCom_2024_Data\Figure 5\peak2peakData_MK_HS.mat')
        dat_P2P = p2p_MK_HS;
    else
        subject = 'Opal';
        load('D:\Figures\VOP_NatCom_2024_Data\Figure 5\peak2peakData_MK_OP.mat')
        mapping_array1_m1_Opal = [98 110 122 101 113 125 74	86;...
            100	112	124	103	115	127	76	88;...
            102	114	126	105	117	66	78	90;...
            104	116	128	107	119	68	80	92;...
            106	118	97	109	121	70	82	94;...
            108	120	99	111	123	72	84	96];
        chnM1 = reshape(mapping_array1_m1_Opal,1,[]);
        dat_P2P = p2p_MK_OP;
    end

    for n = 1:size(dat_P2P,1)
        temp = dat_P2P{n,1}(1:length(chnM1),:);

        % Determine the mean for each channel
        datBarP2P(:,n) = mean(temp,2,'omitnan');
    end

    % Create the box plot traces
    LABEL = ["10 Hz", "50 Hz","80 Hz","100 Hz"];
    F2(monk) = figure;
    F2(monk).Name = ['peak2peakBar_' subject];
    boxplot(datBarP2P,'Labels',LABEL)
    title('All Responses'),ylabel('uV')

    % Run the stats for the animals
alpha = 0.0001;

% Compare 50Hz to all other conditions (two-tail)
idxCond = [1 3 4];
for n = 1:3
    [ciBoot{monk,n,1},rnBoot(n,1)] = bootstrapCompMeans(datBarP2P(:,2),...
        datBarP2P(:,idxCond(n)),10000,alpha,3);
    [ciBoot{monk,n,2},rnBoot(n,2)] = bootstrapCompMeans(datBarP2P(:,2),...
        datBarP2P(:,idxCond(n)),10000,alpha,3,'r');
end

[ciBoot{monk,n,3},rnBoot(n,3)] = bootstrapCompMeans(datBarP2P(:,3),...
    datBarP2P(:,4),10000,alpha,3,'r');
end
