%% This code creates the panels for figure 7g, figures are all scaled to
% largest percent increase.

clear all
close all
clc

PATHNAME = 'D:\Figures\VOP_NatCom_2024\Figure 7';
subj_list = {'S1','S2','S3','S4'};
for s=1:length(subj_list)
    load(fullfile(PATHNAME,strcat([subj_list{s},'.mat'])));
end
S1.hand.aoc80=[];
S1.arm.aoc80=[];
S1.arm.aoc50=[];
S1.arm.aoc100=[];
S1.arm.aocM1=[];
S1.wrist.aoc80=[];
S2.hand.aoc100=[];
S2.wrist.aoc100=[];
S2.arm.aoc80=[];
S2.arm.aoc50=[];
S2.arm.aoc100=[];
S2.arm.aocM1=[];
s=struct;
s.S1=S1;
s.S2=S2;
s.S3=S3;
s.S4=S4;

%% Hand
hand = struct;
subj_list={'S1','S2','S3','S4'};
cc=[0 0 1; 0 1 0; 1 0 0; 1 1 0];
conds = {'aocM1','aoc50','aoc80','aoc100'};

for i=1:length(subj_list)
    data=s.(subj_list{i}).hand;

    % Normalize each data point to a percent of the baseline.
    for c=1:length(conds)
        hand=(data.(conds{c})-nanmean(data.(conds{1})))/nanmean(data.(conds{1}))*100;
        h=nanmean(hand);
        hand_diff(i,c) =h;
    end
end

% Plot
figure
for i=1:length(subj_list)
    scatter([1:4]+0.01*i,hand_diff(i,:),80,'filled','o'), hold on
end
hold on, yline(0)
title('hand')

%% Wrist
wrist = struct;
subj_list={'S1','S2','S3','S4'};
cc=[0 0 1; 0 1 0; 1 0 0; 1 1 0];
conds = {'aocM1','aoc50','aoc80','aoc100'};

for i=1:length(subj_list)
    data=s.(subj_list{i}).wrist;

    % Normalize each data point to a percent of the baseline.
    for c=1:length(conds)
        wrist=(data.(conds{c})-nanmean(data.(conds{1})))/nanmean(data.(conds{1}))*100;
        h=nanmean(wrist);
        wrist_diff(i,c) =h;
    end
end

% Plot
figure
for i=1:length(subj_list)
    scatter([1:4]+0.01*i,wrist_diff(i,:),80,'filled','o'), hold on
end
hold on, yline(0)
title('wrist')

%% Arm
arm = struct;
subj_list={'S1','S2','S3','S4'};
cc=[0 0 1; 0 1 0; 1 0 0; 1 1 0];
conds = {'aocM1','aoc50','aoc80','aoc100'};

for i=1:length(subj_list)
    data=s.(subj_list{i}).arm;

    % Normalize each data point to a percent of the baseline.
    for c=1:length(conds)
        arm=(data.(conds{c})-nanmean(data.(conds{1})))/nanmean(data.(conds{1}))*100;
        h=nanmean(arm);
        arm_diff(i,c) =h;
    end
end

% Plot
figure
for i=1:length(subj_list)
    scatter([1:4]+0.01*i,arm_diff(i,:),80,'filled','o'), hold on
end
hold on, yline(0)
title('arm')