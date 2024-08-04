%% Plot the DCS figures
contactAll = {'APB','BI','DEL','TRI','FLEX','EXT'}; % Best are APB, FLEX, EXT
contactAll = {'APB'};
session = '20220913';
PATHFILE = 'P:\projects\human\VOP STIM\DATA\Intra_Op_DBS\XLTEK';
time = linspace(0,100,600);
xLimPlt = [20 75];
dt = mean(diff(time));

stimCond = {'No stim VIM before','Center, -2.5mm, 50 Hz, 4mA, 100uS',...
    'Anterior, -2.5mm, 50 Hz, 4mA, 100uS','Center, 0.5mm, 50 Hz, 4mA, 100uS',...
    'Anterior, 0.5mm, 50 Hz, 4mA, 100uS','No stim VIM after'};

allCond = {'3mA'};
for n = 1:length(contactAll)
    contact = contactAll(n)
    
    %% Collect the condition data indices
    % Baseline before
    setM1 = 51:79;
    stimPol{1,1} = setM1;
    
    % Center below bottom tract
    setC_n2p5B = 80:110;
    stimPol{2,1} = setC_n2p5B;
    
    % Anterior below bottom tract
    setA_n2p5B = 111:141;
    stimPol{3,1} = setA_n2p5B;
    
    % Center bottom tract
    setC_n0p5B = 142:172;
    stimPol{4,1} = setC_n0p5B;
    
    % Anterior bottom tract
    setA_n0p5B = 173:203;
    stimPol{5,1} = setA_n0p5B;
    
    % Baseline after
    setM1a = 204:233;
    stimPol{6,1} = setM1a;    
    
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
    subset = 125:400;
    for k = 1:size(stimPol,1)
        aoc{k} = nan(maxArea,1);
        aoc{k}(1:lenVal(k),1) = trapz(time(subset),abs(dat{k}(subset,:)));
    end
    
    %% Plot the data
    F(n,1) = figure('Position',[100 100 800 850]); % Waveform figure
    F(n,1).Name = sprintf('%s_traces_%s_m1',session,contact{1});
    F(n,2) = figure('Position',[100 100 800 850]); % Time trace figure
    F(n,2).Name = sprintf('%s_polarity_time_%s_m1',session,contact{1});

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
    myboxplot(aoc,'box')
    F(n,3) = gcf;
    F(n,3).Name = sprintf('%s_auc_%s_m1',session,contact{1});
    title(sprintf('Area under the curve: %s',contact{1}))
    xlabel('Conditions'), ylabel('Area under the curve')
     set(gca,'XTick',1:numel(stimCond))
    set(gca,'XTickLabels',stimCond,'XTickLabelRotation',30)
end