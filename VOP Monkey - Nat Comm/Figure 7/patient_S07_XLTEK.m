% Patient 20240215
clearvars -except AllData
% close all

%% Arm Analysis
% Plot the DCS figures
area = 'arm';% select if u want to run for face or arm

if strcmp(area,'arm')
contactAll = {'L DEL','L BI','L TRI','L FLEX','L EXT','L APB'}; % Best are APB, FLEX, EXT\
contactAll = {'APB','FLEX'}; % Best are APB, FLEX, EXT
testType = {'DCS C'}; %A is face, B is arm\
subset = 70:600; %Here change according to the plot in line 137
time = linspace(0,100,600); %10 ms/div %The middle number is the ms/div*10
dt = mean(diff(time));
xLimPlt = [1 100];
end 

if strcmp(area,'face')
contactAll = {'L TEMP','L MASS','L ORIS','L MYLO','B CRICO'}; % Best are APB, FLEX, EXT
% contactAll = {'L ORIS'}; % Best are APB, FLEX, EXT
testType = {'DCS A'}; %A is face, C is arm
subset = 70:600;
time = linspace(0,50,600); %Change this
dt = mean(diff(time));
xLimPlt = [1 50];
end

session = '20240215';
PATHFILE = 'P:\projects\human\VOP STIM\DATA\Intra_Op_DBS\XLTEK'; %'C:\Users\emg27\Dropbox\PostDoc\ElviraMarco\IntraOperativeForms\XLTEK\';


stimCond = {'BL',...
    '50Hz VIM',...
    '100Hz VIM',...    
    '50Hz VOP',...
    '100Hz VOP',...
    '50Hz VOA',...
    '100Hz VOA',...    
    'BL after'};

%% Collect the condition data indices, 20uV/div
    
    if strcmp(area,'face')
    % Baseline FACE 
    setM1 = 154:184;
    stimPol{1,1} = setM1; 
    % VOA at 50
    set2 = 323:353;
    stimPol{2,1} = set2;   
    % VOA at 100
    set3 = 354:384;
    stimPol{3,1} = set3;   
    % VOP at 50
    set4 = 247:277;
    stimPol{4,1} = set4;   
    % VOP at 100
    set5 = 289:322;
    stimPol{5,1} = set5;
    % VIM at 50
    set6 = 185:215;
    stimPol{6,1} = set6;
    % VIM at 100
    set6 = 216:246;
    stimPol{7,1} = set6;
    % VIM at 130
    set6 = 385:415;
    stimPol{8,1} = set6;
    % BL after
    set7 = 416:446;
    stimPol{9,1} = set7;
    end

    if strcmp(area,'arm')
    % Baseline ARM
    setM1 = 137:176;
    stimPol{1,1} = setM1; 
    % VIM at 50
    set2 = 218:258;
    stimPol{2,1} = set2;   
    % VIM at 100
    set3 = 177:217;
    stimPol{3,1} = set3;   
    % VOP at 50
    set4 = 300:340;
    stimPol{4,1} = set4;   
    % VOP at 100
    set5 = 259:299;
    stimPol{5,1} = set5;
    % VOA at 50
    set6 = 382:422;
    stimPol{6,1} = set6;
    % VOA at 100
    set6 =341:381;
    stimPol{7,1} = set6;
    % BL after
    set7 = 423:447;
    stimPol{8,1} = set7;
    end
    


allCond = {'10mA'};
%saveDir = uigetdir;
for n = 1:length(contactAll)
    contact = contactAll(n)
    
    
% % Here select which conditions you want to plot %at the beginning keep
% commented.
%      conds_todo = [1,2,3];
%      if n==1
%      for i=1:length(conds_todo)
%       stimPol2{i,1}=stimPol{conds_todo(i)};  
%       stimCond2{i}=stimCond{conds_todo(i)};  
%      end
%      stimPol=stimPol2;
%      stimCond=stimCond2;
%      end
    
    
    %% Load the data and determine the area under the curve values
    [dat aoc] = deal(cell(size(stimPol)));
    lenVal = nan(size(stimPol));
    
    for k = 1:size(stimPol,1)
        dat{k} = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',stimPol{k,1}),'testType',testType);
        lenVal(k) = size(dat{k},2);
    end
    
    maxArea = max(lenVal);
    % Calculate the area under the curve
    for k = 1:size(stimPol,1)

%         if k==2 This is to change time or subset in different conditions!
%             time = 0:200:600;
%             subset=(100,300);
%         else 
%             subset()
%         end

     
        aoc{k} = nan(maxArea,1);
        aoc{k}(1:lenVal(k),1) = trapz(time(subset),abs(dat{k}(subset,:)));
              figure,plot((dat{k}))
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
    myboxplot(aoc,'box')
    title(sprintf('Area under the curve: %s',contact{1}))
    xlabel('Conditions'), ylabel('Area under the curve')
    set(gca,'XTick',1:numel(stimCond))
    set(gca,'XTickLabels',stimCond,'XTickLabelRotation',30)

%      AllData.s20230627.voa.(contactAll{n}(end-1:end))={aoc{1},aoc{2},NaN,NaN,NaN,NaN};

end
%saveFigurePDF(F,saveDir);