%% Plot the DCS figures
contactAll = {'APB','EXT'}; 
session = '20221108';
PATHFILE = 'P:\projects\human\VOP STIM\DATA\Intra_Op_DBS\XLTEK\'
time = linspace(0,100,600);
dt = mean(diff(time));
xLimPlt = [1 100];
stimCond = {'BL','Anterior','Posterior'};
    
    %% Collect the condition data indices
    % Baseline 
    setM1 = 181:210;
    stimPol{1,1} = setM1;
     % Anterior
    setC_n0p5B = 271:300;
    stimPol{2,1} = setC_n0p5B;
      
    % Posterior
    setA_n2p5B = 241:270;
    stimPol{3,1} = setA_n2p5B;
allCond = {'13mA'};
%saveDir = uigetdir;
for n = 1:length(contactAll)
    contact = contactAll(n)
    %% Load the data and determine the area under the curve values
    [dat aoc] = deal(cell(size(stimPol)));
    lenVal = nan(size(stimPol));
    
    for k = 1:size(stimPol,1)
        dat{k} = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',stimPol{k,1}));
        lenVal(k) = size(dat{k},2);
    end
    
    maxArea = max(lenVal);
    % Calculate the area under the curve
    subset = 100:500;
    for k = 1:size(stimPol,1)
        aoc{k} = nan(maxArea,1);
        aoc{k}(1:lenVal(k),1) = trapz(time(subset),abs(dat{k}(subset,:)));
    end
    
    %% Plot the data
    F(n,1) = figure('Position',[100 100 800 850]); % Waveform figure
    F(n,1).Name = sprintf('%s_traces_%s_m1',session,contact{1});
    F(n,2) = figure('Position',[100 100 800 850]); % Time trace figure
    F(n,2).Name = sprintf('%s_polarity_time_%s_m1',session,contact{1});
    F(n,3) = figure('Position',[100 100 800 850]); % DCS figure
    F(n,3).Name = sprintf('%s_auc_%s_m1',session,contact{1});
    for k = 1:size(stimPol,1)
        % Plot the raw waveforms
        figure(F(n,1))
        axeD(k) = subplot(size(stimPol,1),2,2*(k-1)+1); hold on
        plot(time,dat{k})
        plot(time,mean(dat{k},2),'k','LineWidth',1.5)
        xlabel('Time (ms)'), ylabel('Trace')
        title(sprintf('DCS %s \n %s',contact{1},stimCond{k}))
        
        % Plot the rectified waveforms
        figure(F(n,1))
        axeR(k) = subplot(size(stimPol,1),2,2*k); hold on
        plot(time,abs(dat{k}))
        plot(time,mean(abs(dat{k}),2),'k','LineWidth',1.5)
        xlabel('Time (ms)'), ylabel('Trace')
        title(sprintf('Rectified DCS %s \n %s',contact{1},stimCond{k}))
        
        % Plot the histogram of the area under the curve
        figure(F(n,2))
        axeT(1) = subplot(size(stimPol,1),1,k); hold on
        d = dat{k}(subset,:);
        if ~isempty(d)
            plot(0:dt:dt*[length(d(:))-1],d(:))
        end
        title(sprintf('DCS through time %s \n %s',contact{1},stimCond{k}))
        xlabel('Time (ms)'), ylabel('Trace')
    end
    
    % Plot the histogram of the area under the curve
    figure(F(n,3))
    %boxplot([aoc{:}])
    myboxplot(aoc,'box')
    title(sprintf('Area under the curve: %s',contact{1}))
    xlabel('Conditions'), ylabel('Area under the curve')
     set(gca,'XTick',1:numel(stimCond))
    set(gca,'XTickLabels',stimCond,'XTickLabelRotation',30)
        
end