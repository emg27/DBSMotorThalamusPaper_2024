%% Plot the DCS figures
contactAll = {'APB','FLEX','EXT'}; 
session = '20220602';
PATHFILE = 'P:\projects\human\VOP STIM\DATA\Intra_Op_DBS\XLTEK'
time = linspace(0,100,600);
dt = mean(diff(time));
xLimPlt = [20 75];
allCond = {'13mA','14mA','15mA'};

for n = 1:length(contactAll)
    contact = contactAll(n)
    for k = 1:length(allCond)
        cond = allCond{k};
        switch cond
            case '13mA'
                % Contact 3: VOP DBS Stim Off 13mA (166-190)
                setM1 = 166:190;
                % Contact 3, VOP -1 +8 50Hz 13mA 190-220
                set50 = 191:220;% 
                % COntact 3, VOP -1 +8 80Hz, 13mA 290-310
                set80 = 281:310;% 
            case '14mA'
                % Contact 3: VOP DBS Stim Off 14mA (153-165) 
                setM1 = [153:165];
                stimPol{5,1} = [153:165];
                stimPol{6,1} = [571:600];
                % Contact 3, VOP -1 +8 50Hz, 14mA 223-250
                set50 = 221:250;
                stimPol{1,1} = 223:250;                
                % COntact 3, VOP -1 +8 80Hz, 14mA 311:340
                set80 = 311:340;
                % COntact 3, VOP +1 -[234] 50Hz, 14mA 393-420
                stimPol{2,1} = 391:420;
                % COntact 3, VOP +1 -[567] 50Hz, 14mA 457-480
                stimPol{3,1} = 451:480;% 
                % COntact 3, VOP +1 -[567] 50Hz, 14mA 521-540
                stimPol{4,1} = 511:540;% 
            case '15mA'
                % Contact 3: VOP DBS Stim Off 15mA (126-144) and after
                % (541-570)
                setM1 = [126:152];%[126:144];% 541:570];
                stimPol{5,2} = [126:152];%[126:144];
                stimPol{6,2} = [541:570];
                % Contact 3, VOP -1 +8 50Hz, 15mA: 251-272
                set50 = 251:280;% 
                stimPol{1,2} = 251:280;% Updated 251:272;
                % COntact 3, VOP -1 +8 80Hz, 15mA 341:370
                set80 = 341:370;                
                % COntact 3, VOP +1 -[234] 50Hz, 15mA 371-387
                stimPol{2,2} = 371:390;% 
                % COntact 3, VOP +1 -[567] 50Hz, 15mA 421-442
                stimPol{3,2} = 421:450;%
                % COntact 3, VOP +1 -8 50Hz, 15mA 481-506
                stimPol{4,2} = 481:510;% 
        end

        % 10 mA M1 only
        dataM1 = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',setM1'));

        % VOP DBS On (3mA 50 Hz) 10 mA M1
        data50 = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',set50'));

        % VOP DBS On (3mA 100 Hz) 10 mA M1
        data80 = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',set80'));

        lM1= size(dataM1,2);
        l50 = size(data50,2);
        l80 = size(data80,2);

        [aocM1, aoc50, aoc80] = deal(nan(1,max([lM1 l50 l80])));

        Fdcs(k,1) = figure('Position',[100 100 800 850]);
        axeD(1) = subplot(3,2,1); hold on
        plot(time,dataM1)
        plot(time,mean(dataM1,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n M1 stim only %s',contact{1},cond))

        axeR(1) = subplot(3,2,2); hold on
        plot(time,abs(dataM1))
        plot(time,mean(abs(dataM1),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n M1 stim only %s',contact{1},cond))

        axeD(2) = subplot(3,2,3); hold on
        plot(time,data50)
        plot(time,mean(data50,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n VOP 50Hz + M1 %s',contact{1},cond))

        axeR(2) = subplot(3,2,4); hold on
        plot(time,abs(data50))
        plot(time,mean(abs(data50),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n VOP 50Hz + M1 %s',contact{1},cond))

        axeD(3) = subplot(3,2,5); hold on
        plot(time,data80)
        plot(time,mean(data80,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n VOP 80Hz + M1 %s',contact{1},cond))
        xlabel('Time (ms)'), ylabel('Trace')

        axeR(3) = subplot(3,2,6); hold on
        plot(time,abs(data80))
        plot(time,mean(abs(data80),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n VOP 80Hz + M1 %s',contact{1},cond))

        matchAxis(axeD)
        setLimits.xLim = xLimPlt;
        matchAxis(axeD,'setLimits',setLimits,'setY',0,'setZ',0,'setR',0,'setTheta',0)
        matchAxis(axeR)
        matchAxis(axeR,'setLimits',setLimits,'setY',0,'setZ',0,'setR',0,'setTheta',0)
        Fdcs(k,1).Name = sprintf('%s_dcs_traces_%s_m1_%s',session,contact{1},cond);

        % Subselect the indices to include in the area under the curve
        subset = 125:400;

        aocM1(1:lM1) = trapz(time(subset),abs(dataM1(subset,:)));
        aoc50(1:l50) = trapz(time(subset),abs(data50(subset,:)));
        aoc80(1:l80) = trapz(time(subset),abs(data80(subset,:)));

        % Plot the histogram of the area under the curve
        myboxplot({aocM1; aoc50; aoc80},'box')
        Fdcs(k,2) = gcf;
        title(sprintf('Area under the curve: %s \n Stim Settings: M1 %s',contact{1},cond))
        xlabel('Conditions'), ylabel('Area under the curve')
        set(gca,'XTick',1:3,'XTickLabels',{'M1 Only','VOP 50 Hz','VOP 80 Hz'})
        Fdcs(k,2).Name = sprintf('%s_dcs_auc_%s_m1_%s',session,contact{1},cond);

        % Plot trace response through time
        Fdcs(k,3) = figure;%('Position',[100 100 800 850]);
        axeT(1) = subplot(3,1,1);
        d = dataM1(subset,:);
        plot(0:dt:dt*[length(d(:))-1],d(:))
        xlabel('Time (ms) [Concate]'), ylabel('Trace')
        title(sprintf('Rect. DCS through time %s \n M1 stim only %s',contact{1},cond))

        axeT(2) = subplot(3,1,2);
        d = data50(subset,:);
        plot(0:dt:dt*[length(d(:))-1],d(:))
        xlabel('Time (ms) [Concate]'), ylabel('Trace')
        title(sprintf('Rect. DCS through time %s \n VOP 50Hz + M1 %s',contact{1},cond))

        axeT(3) = subplot(3,1,3);
        d = data80(subset,:);
        plot(0:dt:dt*[length(d(:))-1],d(:))
        xlabel('Time (ms) [Concate]'), ylabel('Trace')
        title(sprintf('Rect. DCS through time %s \n VOP 80Hz + M1 %s',contact{1},cond))
        
        matchAxis(axeT)
        Fdcs(k,3).Name = sprintf('%s_dcs_time_%s_m1_%s',session,contact{1},cond);
    end
end