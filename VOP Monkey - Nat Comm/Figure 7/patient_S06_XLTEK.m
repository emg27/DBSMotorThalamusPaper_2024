%% Plot the DCS figures
contactAll = {'R DEL','R BI','R TRI','R FLEX','R EXT','R APB'}; 
contactAll = {'R FLEX','R EXT','R APB'}; 
session = '20221108';
PATHFILE = 'P:\projects\human\VOP STIM\DATA\Intra_Op_DBS\XLTEK\'
time = linspace(0,100,600);
dt = mean(diff(time));
xLimPlt = [1 100];
allCond = {'Distance-5.7'};
for n = 1:length(contactAll)
    contact = contactAll(n)
    for k = 1:length(allCond)
        cond = allCond{k};
        switch cond
            case 'Distance-5.7'
                setBL = 181:210;
               
                setVOA = 271:300;
                
                setVOP = 211:240;%
                
                setVIM = 241:270;%
            case 'Distance-3.7'
                setBL = 181:210;
                % Contact 3: VOP DBS Stim Off 14mA (153-165) and after
                % (571-600)
                setVOA =361:390;% 571:600];
              
                setVOP = 301:330;% 
               
                setVIM = 331:360;
        end

        % 10 mA M1 only
        dataBL = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',setBL'));

        % 10 mA M1 only
        dataVOA = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',setVOA'));

        % VOP DBS On (3mA 50 Hz) 10 mA M1
        dataVOP = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',setVOP'));

        % VOP DBS On (3mA 100 Hz) 10 mA M1
        dataVIM = XLTEK.export(session,PATHFILE,'contact',contact,'setNum',...
            compose('#%d',setVIM'));

        lBL= size(dataBL,2);
        lVOA= size(dataVOA,2);
        lVOP = size(dataVOP,2);
        lVIM = size(dataVIM,2);

        [aocBL, aocVOP, aocVIM, aocVOA] = deal(nan(1,max([lBL lVOP lVIM lVOA])));
        %dataM1 = abs(dataM1);
        %data50 = abs(data50);
        %data100 = abs(data100);

        Fdcs(k,1) = figure('Position',[100 100 800 850]);
        axeD(1) = subplot(3,2,1); hold on
        plot(time,dataVOA)
        plot(time,mean(dataVOA,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n VOA %s',contact{1},cond))

        axeR(1) = subplot(3,2,2); hold on
        plot(time,abs(dataVOA))
        plot(time,mean(abs(dataVOA),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n VOA %s',contact{1},cond))

        axeD(2) = subplot(3,2,3); hold on
        plot(time,dataVOP)
        plot(time,mean(dataVOP,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n VOP %s',contact{1},cond))

        axeR(2) = subplot(3,2,4); hold on
        plot(time,abs(dataVOP))
        plot(time,mean(abs(dataVOP),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n VOP  %s',contact{1},cond))

        axeD(3) = subplot(3,2,5); hold on
        plot(time,dataVIM)
        plot(time,mean(dataVIM,2),'k','LineWidth',1.5)
        title(sprintf('DCS %s \n VIM %s',contact{1},cond))
        xlabel('Time (ms)'), ylabel('Trace')

        axeR(3) = subplot(3,2,6); hold on
        plot(time,abs(dataVIM))
        plot(time,mean(abs(dataVIM),2),'k','LineWidth',1.5)
        title(sprintf('Rectified DCS %s \n VIM %s',contact{1},cond))

        matchAxis(axeD)
        setLimits.xLim = xLimPlt;
        matchAxis(axeD,'setLimits',setLimits,'setY',0,'setZ',0,'setR',0,'setTheta',0)
        matchAxis(axeR)
        matchAxis(axeR,'setLimits',setLimits,'setY',0,'setZ',0,'setR',0,'setTheta',0)
        Fdcs(k,1).Name = sprintf('%s_dcs_traces_%s_m1_%s',session,contact{1},cond);

        % Subselect the indices to include in the area under the curve
        subset = 100:500;%1:600;%150:300;

        aocBL(1:lBL) = trapz(time(subset),abs(dataBL(subset,:)));
        aocVOP(1:lVOP) = trapz(time(subset),abs(dataVOP(subset,:)));
        aocVIM(1:lVIM) = trapz(time(subset),abs(dataVIM(subset,:)));
        aocVOA(1:lVOA) = trapz(time(subset),abs(dataVOA(subset,:)));

        aoc{1,1}=aocBL';aoc{2,1}=aocVOP';aoc{3,1}=aocVIM';aoc{4,1}=aocVOA';
%         %PEAK TO PEAK HERE
%         aocM1(1:lM1) = max(dataVOA(subset,:))-min(dataVOA(subset,:));
%         aoc50(1:l50) = max(dataVOP(subset,:))-min(dataVOP(subset,:));
%         aoc100(1:l100) = max(dataVIM(subset,:))-min(dataVIM(subset,:));

        myboxplot(aoc,'box')
        Fdcs(k,2) = gcf;
        title(sprintf('Area under the curve: %s \n Stim Settings: M1 %s',contact{1},cond))
        xlabel('Conditions'), ylabel('Area under the curve')
        set(gca,'XTick',1:numel(aoc))
        set(gca,'XTickLabels',{'M1 Only','VOP 50 Hz','VIM 50 Hz','VOA 50 Hz'})
        Fdcs(k,2).Name = sprintf('%s_dcs_auc_%s_m1_%s',session,contact{1},cond);
    end
end