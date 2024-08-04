contactAll = {'APB', 'FLEX','EXT'}; % Best are APB, FLEX, EXT
session = '20220317';
%PATHFILE = 'D:\Dropbox\PostDoc\ElviraMarco\IntraOperativeForms\XLTEK';
PATHFILE = 'P:\projects\human\VOP STIM\DATA\Intra_Op_DBS\XLTEK'

allCond = {'9mA'};%{'10mA','9mA','8mA'};
time = linspace(0,100,600);
dt = mean(diff(time));
xLimPlt = [10 60];
for n = 1:length(contactAll)
    contact = contactAll(n)
    for k = 1:length(allCond)
        cond = allCond{k};
        switch cond
            case '10mA'
                setM1 = {'#12','#13','#14','#15','#16','#17','#18','#19','#20','#21'};
                set50 = {'#35','#36','#37','#38','#39','#40','#41','#42','#43','#44'};
                set100 = {'#59','#60','#61','#62','#63','#64','#65','#66','#67'};
            case '9mA'
                setM1 = {'#22','#23','#24'};
                set50 = {'#46','#47','#48','#49','#50','#51'};
                set100 = {'#52','#53','#54','#55','#56','#57','#58'};
            case '8mA'
                setM1 = {'#25','#26','#27','#28','#29'};
                set50 = {'#68','#69','#70','#71'};
                set100 = {'#72','#73','#74','#75','#76'};
            case '7mA'
                setM1 = {'#30','#31'};
        end
        lM1= length(setM1);
        l50 = length(set50);
        l100 = length(set100);

        [aocM1, aoc50, aoc100] = deal(nan(1,max([lM1 l50 l100])));

        % 10 mA M1 only
        dataM1 = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',setM1);

        % VOP DBS On (3mA 50 Hz) 10 mA M1
        data50 = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',set50);

        % VOP DBS On (3mA 100 Hz) 10 mA M1
        data100 = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',set100);

        F(k,n) = figure('Position',[100 100 800 850]);
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
        plot(time,data100)
        plot(time,mean(data100,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n VOP 100Hz + M1 %s',contact{1},cond))
        xlabel('Index'), ylabel('Trace')

        axeR(3) = subplot(3,2,6); hold on
        plot(time,abs(data100))
        plot(time,mean(abs(data100),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n VOP 100Hz + M1 %s',contact{1},cond))

        matchAxis(axeD)
        setLimits.xLim = xLimPlt;
        matchAxis(axeD,'setLimits',setLimits,'setY',0,'setZ',0,'setR',0,'setTheta',0)
        matchAxis(axeR)
        matchAxis(axeR,'setLimits',setLimits,'setY',0,'setZ',0,'setR',0,'setTheta',0)
        F(k,n).Name = sprintf('%s_dcs_traces_%s_m1_%s',session,contact{1},cond);

        % Subselect the indices to include in the area under the curve
        subset = 125:400;%150:300;

        aocM1(1:lM1) = trapz(time(subset),abs(dataM1(subset,:)));
        aoc50(1:l50) = trapz(time(subset),abs(data50(subset,:)));
        aoc100(1:l100) = trapz(time(subset),abs(data100(subset,:)));

        % Plot the histogram of the area under the curve
        myboxplot({aocM1; aoc50; aoc100},'box')
        Fdcs(k,n) = gcf;
        title(sprintf('Area under the curve: %s \n Stim Settings: M1 %s',contact{1},cond))
        xlabel('Conditions'), ylabel('Area under the curve')
        set(gca,'XTick',1:3,'XTickLabels',{'M1 Only','VOP 50 Hz','VOP 100 Hz'})
        Fdcs(k,n).Name = sprintf('%s_dcs_auc_%s_m1_%s',session,contact{1},cond);
    end
end