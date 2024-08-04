%% Plot the DCS figures
contactAll = {'APB','BI','DEL','TRI','FLEX','EXT'}; % Best are APB, FLEX, EXT
% contactAll = {'TRI'}; 
session = '20220630';
PATHFILE = 'P:\projects\human\VOP STIM\DATA\Intra_Op_DBS\XLTEK';
time = linspace(0,100,600);
dt = mean(diff(time));
xLimPlt = [10 75];
allCond = {'7mA','11mA'};
for n = 1:length(contactAll)
    contact = contactAll(n)
    for k = 1%:length(allCond)
        cond = allCond{k};
         
        switch cond
            case '7mA'
                % Baseline Before
                setM1 = 192:220;
                stimPol{1,1} = setM1;
                % 7 mA M1, 50 Hz DBS (100PW, 3mA)
                set50 = 252:280;
                stimPol{2,1} = set50;
                % 7 mA M1, 80 Hz DBS (100PW, 3mA)
                set80 = 281:310;
                % 7 mA M1, 100 Hz DBS (100PW, 3mA)
                set100 = 311:340;
                % 7 mA M1, 50 Hz DBS (100PW, 3mA,-1/+2,3,4)
                stimPol{3,1} = 441:471;
                % 7 mA M1, 50 Hz DBS (100PW, 3mA,+1/-2,3,4)
                stimPol{4,1} = 472:500;
                % 7 mA M1, 50 Hz DBS (100PW, 3mA,-8/+5,6,7)
                stimPol{5,1} = 501:530;
                % 7 mA M1, 50 Hz DBS (100PW, 3mA,+8/-5,6,7)
                stimPol{6,1} = 531:560;
                % Baseline After
                setM1aft = 561:590;
                stimPol{7,1} = setM1aft;
            case '11mA'
                % Baseline Before
                setM1 = 221:249;
                % 11 mA M1, 50 Hz DBS (100PW, 3mA)
                set50 = 341:370;
                % 11 mA M1, 80 Hz DBS (100PW, 3mA)
                set80 = 371:410;
                % 11 mA M1, 100 Hz DBS (100PW, 3mA)
                set100 = 411:440;
                % Baseline After
                setM1aft = 591:620;
        end

        % M1 only
        dataM1 = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',setM1'));

        % VOP DBS On (3mA 50 Hz)
        data50 = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',set50'));

        % VOP DBS On (3mA 80 Hz)
        data80 = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',set80'));

        % VOP DBS On (3mA 100 Hz)
        data100 = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',set100'));

        % M1 only after
        dataM1aft = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',setM1aft'));

        lM1= size(dataM1,2);
        lM1a= size(dataM1aft,2);
        l50 = size(data50,2);
        l80 = size(data80,2);
        l100 = size(data100,2);

        [aocM1,aoc50,aoc80,aoc100,aocM1aft] = deal(nan(1,max([lM1 l50 l80 l100 lM1a])));

        % Plot the raw traces
        Fdcs(k,1) = figure('Position',[100 100 800 850]);
        axeD(1) = subplot(5,2,1); hold on
        plot(time,dataM1)
        plot(time,mean(dataM1,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n M1 stim only %s',contact{1},cond))

        axeR(1) = subplot(5,2,2); hold on
        plot(time,abs(dataM1))
        plot(time,mean(abs(dataM1),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n M1 stim only %s',contact{1},cond))

        axeD(2) = subplot(5,2,3); hold on
        plot(time,data50)
        plot(time,mean(data50,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n VOP 50Hz + M1 %s',contact{1},cond))

        axeR(2) = subplot(5,2,4); hold on
        plot(time,abs(data50))
        plot(time,mean(abs(data50),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n VOP 50Hz + M1 %s',contact{1},cond))

        axeD(3) = subplot(5,2,5); hold on
        plot(time,data80)
        plot(time,mean(data80,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n VOP 80Hz + M1 %s',contact{1},cond))
        xlabel('Time (ms)'), ylabel('Trace')

        axeR(3) = subplot(5,2,6); hold on
        plot(time,abs(data80))
        plot(time,mean(abs(data80),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n VOP 80Hz + M1 %s',contact{1},cond))

        axeD(4) = subplot(5,2,7); hold on
        plot(time,data100)
        plot(time,mean(data100,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n VOP 100Hz + M1 %s',contact{1},cond))
        xlabel('Time (ms)'), ylabel('Trace')

        axeR(4) = subplot(5,2,8); hold on
        plot(time,abs(data100))
        plot(time,mean(abs(data100),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n VOP 100Hz + M1 %s',contact{1},cond))

        axeD(5) = subplot(5,2,9); hold on
        plot(time,dataM1aft)
        plot(time,mean(dataM1aft,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n M1 stim only (after) %s',contact{1},cond))

        axeR(5) = subplot(5,2,10); hold on
        plot(time,abs(dataM1aft))
        plot(time,mean(abs(dataM1aft),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n M1 stim only (after) %s',contact{1},cond))

        setLimits.xLim = xLimPlt;
        matchAxis(axeD)
        matchAxis(axeD,'setLimits',setLimits,'setY',0,'setZ',0,'setR',0,'setTheta',0)
        matchAxis(axeR)
        matchAxis(axeR,'setLimits',setLimits,'setY',0,'setZ',0,'setR',0,'setTheta',0)
        Fdcs(k,1).Name = sprintf('%s_dcs_traces_%s_m1_%s',session,contact{1},cond);

        % Subselect the indices to include in the area under the curve
        subset = 125:400;%1:600;%150:300;
        
        aocM1(1:lM1) = trapz(time(subset),abs(dataM1(subset,:)));
        aoc50(1:l50) = trapz(time(subset),abs(data50(subset,:)));
        aoc80(1:l80) = trapz(time(subset),abs(data80(subset,:)));
        aoc100(1:l100) = trapz(time(subset),abs(data100(subset,:)));
        aocM1aft(1:lM1a) = trapz(time(subset),abs(dataM1aft(subset,:)));

        %PEAK TO PEAK HERE
        % aocM1(1:lM1) = max(dataM1(subset,:))-min(dataM1(subset,:));
        % aoc50(1:l50) = max(data50(subset,:))-min(data50(subset,:));
        % aoc80(1:l80) = max(data80(subset,:))-min(data80(subset,:));
        % aoc100(1:l100) = max(data100(subset,:))-min(data100(subset,:));
        % aocM1aft(1:lM1a) = max(dataM1aft(subset,:))-min(dataM1aft(subset,:));

        % Plot the histogram of the area under the curve
        myboxplot({aocM1; aoc50; aoc80; aoc100; aocM1aft},'box')
        Fdcs(k,2) = gcf;
        title(sprintf('Area under the curve: %s \n Stim Settings: M1 %s',contact{1},cond))
        xlabel('Conditions'), ylabel('Area under the curve')
        set(gca,'XTick',1:5,'XTickLabels',{'M1 Only','VOP 50 Hz','VOP 80 Hz','VOP 100 Hz','M1 only after'})
        Fdcs(k,2).Name = sprintf('%s_dcs_auc_%s_m1_%s',session,contact{1},cond);

        % Plot trace response through time
        Fdcs(k,3) = figure('Position',[1 49 1536 740.8000]);%('Position',[100 100 800 850]);
        axeT(1) = subplot(5,1,1);
        d = dataM1(subset,:);
        plot(0:dt:dt*[length(d(:))-1],d(:))
        %plot(abs(d(:)))
        xlabel('Time (ms) [Concate]'), ylabel('Trace')
        title(sprintf('Rect. DCS through time %s \n M1 stim only %s',contact{1},cond))

        axeT(2) = subplot(5,1,2);
        d = data50(subset,:);
        plot(0:dt:dt*[length(d(:))-1],d(:))
        xlabel('Time (ms) [Concate]'), ylabel('Trace')
        title(sprintf('Rect. DCS through time %s \n VOP 50Hz + M1 %s',contact{1},cond))

        axeT(3) = subplot(5,1,3);
        d = data80(subset,:);
        plot(0:dt:dt*[length(d(:))-1],d(:))
        xlabel('Time (ms) [Concate]'), ylabel('Trace')
        title(sprintf('Rect. DCS through time %s \n VOP 80Hz + M1 %s',contact{1},cond))
        
        axeT(4) = subplot(5,1,4);
        d = data100(subset,:);
        plot(0:dt:dt*[length(d(:))-1],d(:))
        xlabel('Time (ms) [Concate]'), ylabel('Trace')
        title(sprintf('Rect. DCS through time %s \n VOP 100Hz + M1 %s',contact{1},cond))

        axeT(5) = subplot(5,1,5);
        d = dataM1aft(subset,:);
        plot(0:dt:dt*[length(d(:))-1],d(:))
        xlabel('Time (ms) [Concate]'), ylabel('Trace')
        title(sprintf('Rect. DCS through time %s \n M1 stim only after %s',contact{1},cond))

        matchAxis(axeT)
        Fdcs(k,3).Name = sprintf('%s_dcs_time_%s_m1_%s',session,contact{1},cond);
    end
end