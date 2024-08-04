% This code creates the bar plot for figure 1b, 1d, 6a, and 7b 

clear 
close all

%% MK-SC

% val = [0 1331 1725];
% vll = [56 7638 1270];
% vpl = [388 8 10];

% val = [0 519 334 330];
% vll = [196 1088 100 102];
% vpl = [356 246 0 0];

% 1-Scorpion 2-Subzero 3-Opal 4-Hosu

mk = {'scorpion','subzero','opal','hosu'};

% val = [0 55.9062 57.875 123.484
%        0 42.9062 42.3281 27
%        21.7188 88.3594 50.2969 40.2188
%        0 15.7031 17.5625 36.7344];
% vll = [54.8438 128.25 53.4062 67.5312
%        0 36.1406 79.4688 105.078
%        0 149.844 76.9844 112.75
%        94.8906 203.047 35.3906 85.9375];
% vpl = [72.4844 70.0156 0 0
%        0 65.2812 36.1406 0
%        165.188 220.703 69.7969 123.422
%        178.891 234.484 18.2812 10.8125];

% Paper Values - Identified using DTI Studio
vpl = [72.4844 70.0156 0 0
       165.188 220.703 69.7969 123.422
       178.891 234.484 18.2812 10.8125];

vll = [54.8438 128.25 53.4062 67.5312
       0 149.844 76.9844 112.75
       94.8906 203.047 35.3906 85.9375];

val = [0 55.9062 57.875 123.484
       21.7188 88.3594 50.2969 40.2188
       0 15.7031 17.5625 36.7344];

% Reviewer Values (25-250mm)
% vpl = [36.6875 14.4219 6.03125 12.3594
%        72.4531 109.781 62.4062 43.1562
%        126.656 167.141 5.0625 16.2344];
% 
% vll = [24.1719 81.4688 40.9531 68.6406
%        25.5 112.203 70.6719 88.5312
%        43.875 167.156 44.2812 75.9219];
% 
% val = [0 54.1719 54.5312 72.4375
%        41.3281 97.2969 33.3594 63.7031
%        31.5781 11.375 80.4375 37.1094];

%% Normalize the area of each region by animal, creates Figure 1b
mthal = {vpl,vll,val};


for ii = 1:size(vpl,1) % # animals
    region{ii}(1,:) = vpl(ii,:)/sum(vpl(ii,:));
    region{ii}(2,:) = vll(ii,:)/sum(vll(ii,:));
    region{ii}(3,:) = val(ii,:)/sum(val(ii,:));
end


nucleus = {};
for nuc = 1:3
    for animal = 1:size(vpl)
        nucleus{nuc}(animal,:) = mthal{1,nuc}(animal,:)/sum(mthal{1,nuc}(animal,:));
    end
end


meanvpl = mean(nucleus{1,1},1)/sum(mean(nucleus{1,1},1));
meanvll = mean(nucleus{1,2},1)/sum(mean(nucleus{1,2},1));
meanval = mean(nucleus{1,3},1)/sum(mean(nucleus{1,3},1));

meanregion(1,:) = meanvpl;
meanregion(2,:) = meanvll;
meanregion(3,:) = meanval;

norm = [meanvpl meanvll meanval];
stdError = [];



for thal = 1:3
    for cort = 1:4
        sc = region{1}(thal,cort);
        op = region{2}(thal,cort);
        hs = region{3}(thal,cort);
        
        error = std([sc op hs])/(3^0.5);
        stdError(thal,cort) = error;
    end
end

figure;
b = bar(meanregion,'grouped');

b(1).FaceColor='y';
b(2).FaceColor='b';
b(3).FaceColor='g';
b(4).FaceColor='m';

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(meanregion);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',meanregion,stdError,'k','linestyle','none');
for reg = 1:size(region,2)
    for animal = 1:size(region{1,1},1)
        s = swarmchart(x(:,reg)',nucleus{1,reg}(animal,:),1000,'red','.','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
        
    end
end

title('mean volume')
hold off

%% PMv

x = [1 2 3;1 2 3; 1 2 3; 1 2 3];

% %PMv - VAL, VLL, VPL
% sc_cst = [];
sz_cst = [47 3.6875 5.23438];
op_cst = [6.04688 17.875 0];
hs_cst = [98.9531 80.875 6.21875];
jc_cst = [24.5781 36.3594 24.4531];

% sc_norm = sc_cst/sum(sc_cst);
sz_norm = sz_cst/sum(sz_cst);
op_norm = op_cst/sum(op_cst);
hs_norm = hs_cst/sum(hs_cst);
jc_norm = jc_cst/sum(jc_cst);


all = [hs_norm;sz_norm;op_norm;jc_norm];
cst_mean = mean(all);
error = std(all)/(4^0.5);

figure;
b = bar(cst_mean);

hold on
% Calculate the number of groups and number of bars in each group

% Plot the errorbars
errorbar([1 2 3],cst_mean,error,'k','linestyle','none');
f = swarmchart(x,all,1000,'red','.');
% f.XJitter = 'rand';
% f.XJitterWidth = 0.2;
title('PMv Mthal Distribution')
hold off


%% Tract of normalized HDFT projections from the VLL, creates Figure 1d

xVOP = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4];

opVOP = [0 100.594 54.4744 60.875];
scVOP = [5.10938 157.922 37.6719 75.7656];
hsVOP = [57.9062 152.359 23.2188 51.6719];
jcVOP = [7.625 41.4375 56.1406 59.4375];

normOP = opVOP/sum(opVOP);
normSC = scVOP/sum(scVOP);
normHS = hsVOP/sum(hsVOP);
normJC = jcVOP/sum(jcVOP);

figure;
bar(mean([normOP; normSC; normHS; normJC]))


vopMean = mean([normOP;normSC;normHS;normJC]);

vopError = std([normSC;normHS;normOP;normJC])/(4^0.5);

hold on
% Calculate the number of groups and number of bars in each group

% Plot the errorbars
errorbar([1 2 3 4],vopMean,vopError,'k','linestyle','none');
j = swarmchart(xVOP,[normSC;normHS;normOP;normJC],1000,'red','.');
% j.XJitter = 'rand';
% j.XJitterWidth = 0.1;
title('VLL Projection Volumes')

hold off

%% Human STIM Projections from the VOP/VIM, creates Figure 7b

xVOP = [1 2 3 4; 1 2 3 4; 1 2 3 4];

s04VOP = [968 2029 106 1400];
s01VOP = [1344 2858 0 1068];
s02VOP = [2519 4197 2332 3624];


norm_s04 = s04VOP/sum(s04VOP);
norm_s01 = s01VOP/sum(s01VOP);
norm_s02 = s02VOP/sum(s02VOP);

all = [norm_s04;norm_s01;norm_s02];


figure;
bar(mean(all))


vopMean = mean(all);

vopError = std(all)/(3^0.5);

hold on
% Calculate the number of groups and number of bars in each group

% Plot the errorbars
errorbar([1 2 3 4],vopMean,vopError,'k','linestyle','none');
j = swarmchart(xVOP,all,1000,'red','.');
% j.XJitter = 'rand';
% j.XJitterWidth = 0.1;
title('Human VOP/VIM STIM Projection Volumes')

hold off

%% POST LESION Projections from the VLL, creates Figure 6a
clear
close all

x = [1 2;1 2; 1 2;1 2];

%M1
% sc_cst = [242.016 256.594];
% sz_cst = [463.234 35.4688];
% op_cst = [204.703 84.5156];
% % hs_cst = [186.734 0];
% jc_cst = [323.2455 132.047];

%S1
% sc_cst = [81.5781 5.9375];
% sz_cst = [322.406 110.453];
% op_cst = [133.547 48.6719];
% % hs_cst = [];
% jc_cst = [238.328 65.5];

% %PMd
% sc_cst = [72.1406 1];
% sz_cst = [423.734 282.422];
% op_cst = [278.422 155.391];
% % hs_cst = [];
% jc_cst = [306.562 24.8281];

% %SMA
% sc_cst = [34.5 1];
% sz_cst = [587.656 1];
% op_cst = [226.219 81.2812];
% % hs_cst = [];
% jc_cst = [543.609 124.609];
% 
% %PMv
% sc_cst = [41.0781 1];
% sz_cst = [107.084 76.2969];
% op_cst = [30.7188 14.2812];
% % % hs_cst = [];
% % %jc_cst = [0 0];


% cumalative area
sc_cst = [218.641 228.234];
sz_cst = [485.469 101.375];
op_cst = [424.094 267.688];
% hs_cst = [186.734 0];
jc_cst = [679.047 230.688];

% sc_norm = sc_cst/sum(sc_cst);
% sz_norm = sz_cst/sum(sz_cst);
% op_norm = op_cst/sum(op_cst);
% % hs_norm = hs_cst/sum(hs_cst);
% jc_norm = jc_cst/sum(jc_cst);

sc_norm = sc_cst/max(sc_cst);
sz_norm = sz_cst/max(sz_cst);
op_norm = op_cst/max(op_cst);
% hs_norm = hs_cst/sum(hs_cst);
jc_norm = jc_cst/max(jc_cst);


all = [sc_norm;sz_norm;op_norm;jc_norm];
cst_mean = mean(all);
error = std(all)/(3^0.5);

[ci95, rejectNull] = bootstrapCompMeans(all(:,1),all(:,2),10000,.05)

figure;

b = bar(cst_mean);


hold on
% Calculate the number of groups and number of bars in each group

% Plot the errorbars
errorbar([1 2],cst_mean,error,'k','linestyle','none');
f = swarmchart(x,all,1000,'red','.');
% f.XJitter = 'rand';
% f.XJitterWidth = 0.2;
title('mean volume')
hold off



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









