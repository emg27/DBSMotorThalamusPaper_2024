clear all
close all

pathName = 'D:\Figures\VOP_NatCom_2024_Data\Sup Fig 6\';%'C:\Users\jch83\OneDrive - University of Pittsburgh\Documents\Research - Proprioception - VOP\Data\VOP-Monkey\Figures\Figure 6\JOHNNY\DLC\';
trialsToLoad = [175 181 180 182];

bodypartlst = {'thumb','index','pinky','wrist_m','wrist_l'};

%% Load Data

for t = 1:length(trialsToLoad)
     dlc_filename = sprintf('MK-JC_DLCpost_PAD%03d',trialsToLoad(t));
    load([pathName dlc_filename '.mat'])
    data{t} = dat;
    datacol = 2;
    col = 1;
    for ii = 1:length(bodypartlst)
        
        trial = data{1,t};
        traces{t}(:,col:col+1) = trial(:,datacol:datacol+1);
        datacol = datacol+3;
        col = col+2;
    end
end


%% Average
labels = {'No Stim', '50 Hz', '80Hz','100 Hz'};
for t = 1:size(traces,2)
   for aa = 1:size(traces{1,t},2)
       if rem(aa,2) == 0
            singleTrace = traces{1,t}(:,aa)*-1;
       else
            singleTrace = traces{1,t}(:,aa);
       end

        [pks,locs] = findpeaks(singleTrace,'MinPeakDistance',15,'MinPeakProminence',1);
%         figure;plot(1:length(singleTrace),singleTrace,locs,pks,'o')
        title(labels{t})

        for ii = 1:length(locs)
            pre =locs(ii)-5;
            post = locs(ii)+5;
            if pre > 0 && post < length(singleTrace)
                window = singleTrace(pre:post)';
                allTraces{t,aa}(ii,:) = window-mean(window(1:2));  
            end

        end
   end


    for aa = 1:size(traces{1,t},2)
        cleanedTraces{t,aa} = allTraces{t,aa};
        mean_signal=mean(cleanedTraces{t,aa},2);
        std_signal=std(mean_signal);
        idxToRemove=find(abs(mean_signal)>std_signal*2);
%         cleanedTraces{t,aa}(idxToRemove,:) = [];
        staTrace{t,aa} =mean(cleanedTraces{t,aa},1);
    end
end

%% Peak to Peak

%scale pixels/mm
scale = [1];

 for t = 1:size(traces,2)
   for aa = 1:size(cleanedTraces(t,:),2)
       for ii  = 1:size(cleanedTraces{t,aa})
           p2p = peak2peak(cleanedTraces{t,aa}(ii,:))/scale;

           amps{t,aa}(ii) = p2p;
       end
       amps{t,aa} = rmoutliers(amps{t,aa});
   end
 end

%% Plot Boxplots

close all
for shift = 1:2
    axis = {'X','Y'};
    
    figure;
    count = 1;
    for aa = shift:2:length(bodypartlst)*2

        c = [];
        grp = [];

        for ii = 1:length(trialsToLoad)
            single = amps{ii,aa}';

            c = [c; single];
            grp = [grp; ones(1,length(single))'*ii];
        end
        
        subplot(1,5,count)
        boxplot(c,grp,'Labels',labels,'symbol','+r')
        title([bodypartlst{count} axis{shift}])


        trialMat = [1 2
            1 3];

        alpha = .05;
        offset = 0;
        for comb = 1:size(trialMat,1)
            trial1 = trialMat(comb,1);
            trial2 = trialMat(comb,2);

            [CI,sig]=bootstrapCompMeans(amps{trial1,aa},amps{trial2,aa},10000,alpha/size(trialMat,1));

            if sig
                yt = get(gca, 'YTick');
                %                 axis([xlim   0  ceil(max(yt)*5)])
                annotation('textbox',[.1 .9 .1 .1],'String',"alpha = "+alpha)
                xt = get(gca, 'XTick');
                hold on
                plot(xt([trial1 trial2]), [1 1]*max(yt)*(1.1+offset), '-k',  mean(xt([trial1 trial2])), max(yt)*(1.15+offset), '*k')
                hold off
            end
            offset = offset+0.001;
        end
        count = count+1;
    end
end

%% Percent Variation
% close all
meanp2pX = [];
meanp2py = [];

for shift = 1:2
    axis = {'X','Y'};
    figure;
    count = 1;
    for aa = shift:2:length(bodypartlst)*2
        for ii = 1:length(trialsToLoad)
            if shift ==1
                meanp2pX(ii,count) = median(amps{ii,aa});
            elseif shift==2
                meanp2pY(ii,count) = median(amps{ii,aa});
            end

        end
        if shift ==1
            subplot(1,5,count)
            baseline = meanp2pX(1,count);
            plot((meanp2pX(:,count)-baseline)/baseline*100,'-o')
            title(['Johnny POST ' bodypartlst{count} axis{shift}])
            xticklabels({'ic alone','50Hz','100Hz'})
            ylabel("Percent Variation")

        elseif shift==2
            subplot(1,5,count)
            baseline = meanp2pY(1,count);
            plot((meanp2pY(:,count)-baseline)/baseline*100,'-o')
            title(['Johnny POST ' bodypartlst{count} axis{shift}])
            xticklabels({'ic alone','50Hz','100Hz'})
            ylabel("Percent Variation")

        end

        count = count+1;
    end
    
end

%% PLOT Traces
bodypart = 'all';

savefile = 0;

plotLst = {};
if strcmp(bodypart,'all')
    for bb = 1:length(bodypartlst)
        body = bodypartlst{bb};
        figure;
        parameters = plot_save_dlc(data,body,savefile);
        plotLst{bb} = parameters;
    end
else
    parameters = plot_save_dlc(data,bodypart,savefile);
    plotLst{1} = parameters;
end



%% FUNCTIONS

function parameters = plot_save_dlc(data,bodypart,save)

nostim = 1;

hz50 = 2;
hz80 = 3;
hz100 = 4;


if bodypart == "thumb"
part_x=2;
part_y=3;
elseif bodypart == "index"
part_x=5;
part_y=6;
elseif bodypart == "pinky"
part_x=8;
part_y=9;
elseif bodypart == "wrist_m"
part_x=11;
part_y=12;
elseif bodypart == "wrist_l"
part_x=14;
part_y=15;
end

    nostim_x = data{1,nostim}(100:end,part_x)-mean(data{1,nostim}(100:end,part_x))+10;
    hz80_x = data{1,hz80}(100:end,part_x)-mean(data{1,hz80}(100:end,part_x))+15;
    hz50_x = data{1,hz50}(100:end,part_x)-mean(data{1,hz50}(100:end,part_x))+20;
    hz100_x = data{1,hz100}(100:end,part_x)-mean(data{1,hz100}(100:end,part_x))+25;

    nostim_y = data{1,nostim}(100:end,part_y)-mean(data{1,nostim}(100:end,part_y))+10;
    hz80_y = data{1,hz80}(100:end,part_x)-mean(data{1,hz80}(100:end,part_x))+15;
    hz50_y = data{1,hz50}(100:end,part_y)-mean(data{1,hz50}(100:end,part_y))+20;
    hz100_y = data{1,hz100}(100:end,part_y)-mean(data{1,hz100}(100:end,part_y))+25;

    tiledlayout(1,2)
    time = (1:length(nostim_x))/30;
   
    time50 = (1:length(hz50_x))/30;
    time80 = (1:length(hz80_x))/30;
    time100 = (1:length(hz100_x))/30;

    parameters = {time,time50,time80,time100,nostim_x,hz50_x,hz80_x,hz100_x,nostim_y,hz50_y,hz80_y,hz100_y,bodypart};
    
%     figure;
    ax1 = nexttile;
    hold on 
    title(bodypart + " X")
    ylabel('amplitude (pixels)')
    xlabel('Time (sec)')
    plot(time,nostim_x,'red')
    
    plot(time50,hz50_x,'cyan')
    plot(time80,hz80_x,'green')
    plot(time100,hz100_x, 'black')
    legend('No Stim','50 Hz','80 Hz','100 Hz')
    hold off

    ax2 = nexttile;
    hold on
    title(bodypart + " Y")
    ylabel('amplitude (pixels)')
    xlabel('Time (sec)')
    plot(time,nostim_y,'red')
    
    plot(time50,hz50_y,'cyan')
    plot(time80,hz80_y,'green')
    plot(time100,hz100_y, 'black')
    legend('No Stim','50 Hz','80 Hz','100 Hz')
    hold off

    linkaxes([ax1,ax2],'xy')
   
    
    if save
        saveas(gca,"20211021_DLC_"+"_"+bodypart)
    end
    
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