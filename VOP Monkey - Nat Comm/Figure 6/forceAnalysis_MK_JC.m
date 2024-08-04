clear all
close all

%% INPUT VARIABLES

TrialsToLoad=[169 185 188 187 189];
analog = 2;

pre = 200;
post = 1500;

pathName = 'P:\projects\monkeys\Acute_Proprioception\ANIMALS\Mk-JC M68-22 Johnny Cage\20221005 - Experiment\DATA\mat\';

channels = [129:148];

% Hz = {'IC alone', 'VOP 50 Hz', 'VOP 80 Hz', 'VOP 100 Hz'};
Hz = {'IC PRE', 'IC POST', 'VOP 50 Hz', 'VOP 80 Hz', 'VOP 100 Hz'};
forces = {};

StimThresh=16000;
Fs = 30000;
%% Load Data

for t=1:length(TrialsToLoad)
    ns5_filename = sprintf('NS5_PAD%04d',TrialsToLoad(t));
    load([pathName ns5_filename '.mat']);

    allData{t} = NS5;

    [ftCal] = forceTransform(allData{t});
    allTransducerDat{t} = ftCal;
    forces{t} = allTransducerDat{1,t}(1:3,:);
    torques{t} = allTransducerDat{1,t}(4:6,:);
end

%% GET STIM TIMES
for t=1:length(TrialsToLoad)
    StimThresh=16000;
    Fs = 30000;

    stim=findAnalogStim(allData{t},StimThresh);
    StimTimes{t} = stim;
    vop = stim;
    StimTimes{t}=StimTimes{t}(StimTimes{t}(:,1)==analog,2);
    %     vop = vop(vop(:,1)==3,2);
    %
    %     if length(unique(stim(:,1))) > 1
    %         StimTimes{t} = StimTimes{t}(StimTimes{t}>vop(1));
    %     end
end

%% STA
skip = 47;
stimsToSkip = 12;
stimToInclude = 22;

stimArtifacts = {};
for t=1:length(TrialsToLoad)
    count=1;
    stimCount = 0;
    for aa = 1:skip:length(StimTimes{1,t})
        stimCount = stimCount+1;
        
        if t>1 && stimCount < stimsToSkip
            continue
        end
% 
%         if t==2 && stimCount >= stimToInclude
%             continue
%         end
        if t==1 && stimCount < stimsToSkip
            continue
        end       

        % Truncate the force information to the desired window size.
        stim = StimTimes{1,t}(aa);
        lowerbound = stim-(pre*Fs/1000);
        upperbound = stim+(post*Fs/1000);
        if upperbound < length(forces{1,t}(1,:)) && lowerbound >= 0
            for ii=1:size(forces{1,t},1)
                
                singleForce = forces{1,t}(ii,:);
                singleTorque = torques{1,t}(ii,:);

                window_f = singleForce(lowerbound:upperbound);
                window_t = singleTorque(lowerbound:upperbound);
                window_f = window_f - mean(window_f(1:5000));
                window_t = window_t - mean(window_t(1:5000));
                indForces{t,ii}(count,:) = window_f;
                indtorques{t,ii}(count,:) = window_t;

                if aa+skip-1 < length(StimTimes{1,t})
                    stimArtifacts{t,ii}(:,count) = (StimTimes{1,t}(aa:aa+skip-1)-StimTimes{1,t}(aa)+1)/30000*1000;
                end

            end

            % Calculate the magnitude of the force and torque
            magForce{t}(count,:) = sum([indForces{t,1}(count,:).^2; indForces{t,2}(count,:).^2; indForces{t,3}(count,:).^2],1).^0.5;
            magTorque{t}(count,:) = sum([indtorques{t,1}(count,:).^2; indtorques{t,2}(count,:).^2; indtorques{t,3}(count,:).^2],1).^0.5;

            p2p_f{t}(count,:) = peak2peak(magForce{t}(count,:));
            area = cumtrapz(magForce{t}(count,:));
            auc_f{t}(count,:) = area(end);

            p2p_t{t}(count,:) = peak2peak(magTorque{t}(count,:));
            area = cumtrapz(magTorque{t}(count,:));
            auc_t{t}(count,:) = area(end);


            count=count+1;
        end
    end
    p2p_f{t} = (p2p_f{t});
    auc_f{t} = (auc_f{t});
    p2p_t{t} = (p2p_t{t});
    auc_t{t} = (auc_t{t});

    staForces{t} = mean(magForce{t});
    staTorques{t}= mean(magTorque{t});
end


%% PLOT the force traces and box plots (Figure 6c)
close all
time = [-pre:1/Fs*1000:post];

c = [.7 .7 .7];
cAvg = [0 0 0];

for ii=1:size(staForces,1)
    % Plot the force traces through time.
    figure;
    axes1_f = [];
    axes1_t = [];
    for t=1:length(TrialsToLoad)

        ax1 = subplot(2,length(TrialsToLoad),t);
        hold on
        plot(time,magForce{t},'color',c,'LineWidth',0.8)
        plot(time,staForces{t},'color',cAvg,'LineWidth',1.5)
        
%         plot(stimArtifacts{t,ii}(:,1),ones(length(stimArtifacts{t,ii}(:,1)))*2.5,'*')
        axes1_f(end+1) = ax1;
        title(['Force ' Hz{t}])
        xlim([-200 1500])
        linkaxes(axes1_f,'xy')
        hold off

        ax1 = subplot(2,length(TrialsToLoad),t+length(TrialsToLoad));
        hold on
        plot(time,magTorque{t},'color',c,'LineWidth',0.8)
        plot(time,staTorques{t},'color',cAvg,'LineWidth',1.5)
        
        title(['Torque ' Hz{t}])
        axes1_t(end+1) = ax1;
        xlim([-200 1500])
        linkaxes(axes1_t,'xy')
        hold off
    end

    % Plot the peak to peak box plots
    figure;
    axes2_f = [];
    axes2_t = [];
    for t=1:length(TrialsToLoad)

        ax2 = subplot(2,length(TrialsToLoad),t);
        boxplot(p2p_f{t})
        axes2_f(end+1) = ax2;
        title(['p2p F ' Hz{t}])
        linkaxes(axes2_f,'xy')

        ax2 = subplot(2,length(TrialsToLoad),t+length(TrialsToLoad));
        boxplot(p2p_t{t})
        axes2_t(end+1) = ax2;
        title(['p2p T ' Hz{t}])
        linkaxes(axes2_t,'xy')
    end

    % Plot the area under the curve (AUC) box plots
    figure;
    axes3_f = [];
    axes3_t = [];
    for t=1:length(TrialsToLoad)

        ax3 = subplot(2,length(TrialsToLoad),t);
        boxplot(auc_f{t})
        axes3_f(end+1) = ax3;
        title(['auc F ' Hz{t}])
        linkaxes(axes3_f,'xy')

        ax3 = subplot(2,length(TrialsToLoad),t+length(TrialsToLoad));
        boxplot(auc_t{t})
        axes3_t(end+1) = ax3;
        title(['auc T ' Hz{t}])
        linkaxes(axes3_t,'xy')
    end

end
%% boostrap

[ci95, rejectNull] = bootstrapCompMeans(auc_f{2}, auc_f{5}, 10000,.001)

%% FUNCTIONS

function[ftCal] = forceTransform(dat_ns5, varargin)
% This is the offset that we found for the force transducer
offset = [-24.3805 -12.5061 -11.1104 -27.8208 68.0922 29.1555]/1000;%[0.1149    1.2351   -0.8644   -0.8464    0.7502   -0.8112];%
% This is the calibration matrix for the force transducer that we are using
% from the Batista lab
ftCalMat = [0.06556   0.00618  -0.07958   3.16157  -0.00419  -3.23805;... 'Fx'
    0.05082  -3.84631   0.00356   1.77163  -0.04008   1.88466;... 'Fy'
    5.41187  -0.09469   5.31854  -0.07083   5.36020  -0.04775;...'Fz'
    0.00085  -0.02107   0.07438   0.00967  -0.07298   0.01042;...'Tx'
    -0.09029   0.00093   0.04433  -0.01821   0.04373   0.01831;...'Ty'
    0.00047  -0.04626   0.00115  -0.04386   0.00090  -0.04413]; %'Tz'
chnIdx = 21:-2:11; % There are the analog channels that store the force transducer gauges.

% Find the index of the analog
analogChn = find(contains({dat_ns5.ElectrodesInfo.Label},'analog'));

% Collect the force data
forceTrans = double(dat_ns5.Data(analogChn(chnIdx),:))/1000;

% Calculate the force and torque magnitudes
ftCal = ftCalMat*(forceTrans-offset'*ones(1,size(forceTrans,2)));

% Calculate the magnitude
% forceMag = (sum(ftCal(1:3,:).^2,1)).^0.5;
% forceMag = filterData(forceMag);

% forceMag = sum(ftCal(1:3,:),1);
% torqueMag = (sum(ftCal(4:6,:).^2,1)).^0.5;
% torqueMag = filterData(torqueMag);

for aa = 1:size(ftCal,1)
    ftCal(aa,:) = filterData(ftCal(aa,:));

end

end


function filtered_LFP = filterData(LFP)
fs_chan = 30000;

filterbands_Line=[58,62];
[z, p, k] = butter(2, filterbands_Line/(fs_chan/2), 'stop');
[sos_line,g_line] = zp2sos(z, p, k, 'down', 'two');

filterbands=[6];
[z,p,k] = butter(2, filterbands/(fs_chan/2), 'low');
[sos,g] = zp2sos(z, p, k, 'down', 'two');


filtered_LFP=filtfilt(sos_line,g_line,double(LFP'))';
filtered_LFP=filtfilt(sos,g,double(filtered_LFP'))';

end

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
