function myFishData = EEG_SpikeFinder(DataFolder,varargin)
    %%%% EEG Spike Finder
    % crobens@mit.edu

    %%%% parse parameters %%%%
    p = inputParser;
    
    % time block settings
    p.addParameter('ExcludeInitialMinutes',     5);     % in minutes
    p.addParameter('ExcludeLastMinutes',        1);     % in minutes
    p.addParameter('startTimeBlockInMinutes',   20);    % in minutes
    p.addParameter('endTimeBlockInMinutes',     5);     % in minutes
    
    % peak finder settings
    p.addParameter('SNThresholdNegative',       5);     % 
    p.addParameter('SNThresholdPositive',       1.5);   % 
    p.addParameter('maxFireRate',               10);    % in Hz

    % filter and noise analysis settings
    p.addParameter('downsamplingFactor',        10);    % 
    p.addParameter('alpha',                     0.5);   % 
    p.addParameter('SGorder',                   3);     % 
    p.addParameter('SGlength',                  31);    % 
   
    % seperate normal and epileptic spikes
    p.addParameter('normalEEGsigma',            2.5);   % 
    
    % burst analysis settings
    p.addParameter('burstThresholdBeta',            1);     % in seconds
    p.addParameter('burstThresholdGamma',           3);     % in seconds

    % plotting and debug settings
    p.addParameter('plotting',                  true);  % 
    p.addParameter('plotDebug',                 true);  % 
    p.addParameter('plotInMinutes',             true);  %
    p.addParameter('pauseEachFile',             true);  %
    
    % export file format (eps or png)
    p.addParameter('exportFormat',              'eps'); % 'eps' or 'png'
    
    p.parse(varargin{:});
    
    ExcludeInitialMinutes   = p.Results.ExcludeInitialMinutes;
    ExcludeLastMinutes      = p.Results.ExcludeLastMinutes;
    startTimeBlockInMinutes = p.Results.startTimeBlockInMinutes; 
    endTimeBlockInMinutes   = p.Results.endTimeBlockInMinutes; 
    
    SNThresholdNegative     = p.Results.SNThresholdNegative;
    SNThresholdPositive     = p.Results.SNThresholdPositive;
    maxFireRate             = p.Results.maxFireRate; 
    
    downsamplingFactor      = p.Results.downsamplingFactor;
    alpha                   = p.Results.alpha;
    SGorder                 = p.Results.SGorder;
    SGlength                = p.Results.SGlength;
    
    normalEEGsigma          = p.Results.normalEEGsigma;
    
    burstThresholdBeta      = p.Results.burstThresholdBeta;
    burstThresholdGamma     = p.Results.burstThresholdGamma;
    
    plotting                = p.Results.plotting;
    plotDebug               = p.Results.plotDebug;
    plotInMinutes           = p.Results.plotInMinutes;
    pauseEachFile           = p.Results.pauseEachFile;
    
    exportFormat            = p.Results.exportFormat;
    
   
    
    
    % create empty containers
    gamma_spikeRateInHz           = [];
    gamma_spikeStrength           = [];
    gamma_spikeDuration           = [];
    gamma_spikeAmplitude          = [];
    exportTable             = table;
    exportTableEnd      = table;
    excelTableName          = 'EEG_SpikeFinder_results';

    %make new analysis folder
    formatOut       = 'yyyy-mm-dd_hh-MM';
    timeString      = datestr(now,formatOut);
    analysisFolder  = [timeString '_analysis'];

    mkdir([DataFolder analysisFolder])

    files = dir([DataFolder '*.abf']);
    for fish=1:length(files)
        myFishData = table;
        myFishDataEnd = table;
        % load data
        display(files(fish).name);
        [data, ~, header]       = abfload([DataFolder files(fish).name]);
        % get sampling frequency
        samplingFrequency   = 1/header.si*1e6; % in Hz 

        % downsample data
        dataDownSampled = imresize(data,1/downsamplingFactor,'box');
        timesInSeconds = 1/samplingFrequency*(0:downsamplingFactor:length(data)-1);

        % run SG filter over data
        dataSGfiltered  = sgolayfilt(dataDownSampled,SGorder,SGlength);

        if plotDebug
            figure(881),clf;
            if plotInMinutes
                subplot(2,1,1)
                hold on
                plot(timesInSeconds/60,dataDownSampled)
                plot(timesInSeconds/60,dataSGfiltered,'-','LineWidth',2)
                hold off
                xlim([5,6])
                title('initial trace')
                legend('downsampled data', 'SG filtered data')
                box on
                xlabel('time (min)')
                ylabel('mV')
                set(gca, 'FontName', 'Arial')
                set(gca,'FontSize', 12);

                subplot(2,1,2)
                hold on
                plot(timesInSeconds/60,dataDownSampled)
                plot(timesInSeconds/60,dataSGfiltered,'-','LineWidth',2)
                hold off
                xlim([timesInSeconds(end)/60-5,timesInSeconds(end)/60-4])
                title('end trace')
                legend('downsampled data', 'SG filtered data')
                box on
                xlabel('time (min)')
                ylabel('mV')
                set(gca, 'FontName', 'Arial')
                set(gca,'FontSize', 12);
            else
                subplot(2,1,1)
                hold on
                plot(timesInSeconds,dataDownSampled)
                plot(timesInSeconds,dataSGfiltered,'-','LineWidth',2)
                hold off
                xlim([5*60,6*60])
                title('initial trace')
                legend('downsampled data', 'SG filtered data')
                box on
                xlabel('time (s)')
                ylabel('mV')
                set(gca, 'FontName', 'Arial')
                set(gca,'FontSize', 12);

                subplot(2,1,2)
                hold on
                plot(timesInSeconds,dataDownSampled)
                plot(timesInSeconds,dataSGfiltered,'-','LineWidth',2)
                hold off
                xlim([timesInSeconds(end)-5*60,timesInSeconds(end)-4*60])
                title('end trace')
                legend('downsampled data', 'SG filtered data')
                box on
                xlabel('time (s)')
                ylabel('mV')
                set(gca, 'FontName', 'Arial')
                set(gca,'FontSize', 12);
                
            end
            
        end

        if plotting   
            figure(882),clf;
            imposedYlim = [-0.4,0.1];
            hold on
            plot(timesInSeconds/60,dataSGfiltered,'-','LineWidth',1)
            plot([ExcludeInitialMinutes,ExcludeInitialMinutes],imposedYlim,'k--','LineWidth',2);
            plot([ExcludeInitialMinutes+startTimeBlockInMinutes,ExcludeInitialMinutes+startTimeBlockInMinutes],imposedYlim,'k--','LineWidth',2);
            plot([timesInSeconds(end)/60-ExcludeLastMinutes,timesInSeconds(end)/60-ExcludeLastMinutes],imposedYlim,'k-.','LineWidth',2);
            plot([timesInSeconds(end)/60-ExcludeLastMinutes-endTimeBlockInMinutes,timesInSeconds(end)/60-ExcludeLastMinutes-endTimeBlockInMinutes],imposedYlim,'k-.','LineWidth',2);
            hold off
            xlim([0.5,timesInSeconds(end)/60])
            ylim(imposedYlim)
            title(['full trace: ' files(fish).name])
            legend('SG filtered data')
            box on
            xlabel('time (min)')
            ylabel('V')
            set(gca, 'FontName', 'Arial')
            set(gca,'FontSize', 12);
            if strcmp(exportFormat,'eps')
                if ispc
                    print('-depsc', '-tiff', '-r300', '-painters', [DataFolder analysisFolder '\WholeTraceOverview_ID_' num2str(fish) '_' files(fish).name '.eps']);
                elseif ismac
                    print('-depsc', '-tiff', '-r300', '-painters', [DataFolder analysisFolder '/WholeTraceOverview_ID_' num2str(fish) '_' files(fish).name '.eps']);
                end
            else
                if ispc
                    print('-dpng', '-r300', '-painters', [DataFolder analysisFolder '\WholeTraceOverview_ID_' num2str(fish) '_' files(fish).name '.png']);
                elseif ismac
                    print('-dpng', '-r300', '-painters', [DataFolder analysisFolder '/WholeTraceOverview_ID_' num2str(fish) '_' files(fish).name '.png']);
                end
            end
            figure(831),clf;
            figure(889),clf;
            figure(832),clf;
            figure(690),clf;
            figure(819),clf;
            figure(77112),clf;

        end


        for dataOrControl = 1:2

            if dataOrControl == 1
                dataMask                = timesInSeconds/60>ExcludeInitialMinutes & ...
                                          timesInSeconds/60<startTimeBlockInMinutes+ExcludeInitialMinutes;
                dataToAnalyze           = dataSGfiltered(dataMask);
                timesInSecondsToAnalyze = timesInSeconds(dataMask);
            elseif dataOrControl == 2
                dataMask                = timesInSeconds/60>timesInSeconds(end)/60-endTimeBlockInMinutes-ExcludeLastMinutes & ...
                                          timesInSeconds/60<timesInSeconds(end)/60-ExcludeLastMinutes;
                dataToAnalyze           = dataSGfiltered(dataMask);
                timesInSecondsToAnalyze = timesInSeconds(dataMask);
            end

            % get noise level
            figure(771),clf;
            dataHistogram   = histogram(dataToAnalyze,'BinWidth',1e-3);
            histValues      = dataHistogram.Values;
            histBinCenters  = dataHistogram.BinEdges(1:end-1)+mean(diff(dataHistogram.BinEdges));

            %fit central noise with gaussian to autothreshold signal and noise
            [~,maxIdx] = max(histValues);
            maxHistValue = mean(maxk(histValues,5));
            dataMask = find(histBinCenters<histBinCenters(maxIdx) & histValues>=0.75*maxHistValue,1,'first');
            leftEdge = histBinCenters(dataMask);
            dataMask = find(histBinCenters>histBinCenters(maxIdx) & histValues>=0.75*maxHistValue,1,'last');
            rightEdge= histBinCenters(dataMask);
            
            
            
            if isempty(rightEdge) || isempty(leftEdge)
                fitMask(maxIdx-5:maxIdx+5) = true;
                pGuess      = [maxHistValue,histBinCenters(maxIdx),5];  % Inital (guess) parameters
                leftEdge    = histBinCenters(maxIdx-5);
                rightEdge   = histBinCenters(maxIdx+5);
            else
                fitMask         = histBinCenters>=leftEdge & histBinCenters<=rightEdge;
                pGuess  = [maxHistValue,histBinCenters(maxIdx),range([leftEdge,rightEdge]/2)];  % Inital (guess) parameters
            end
            
            
            if sum(fitMask)<5 
                fitMask(maxIdx-5:maxIdx+5) = true;
            end
            
            fitfun          = @(p,x) p(1).*exp(-1/2*((x-p(2))/p(3)).^2);

            

            lb      = [0,-inf,0];  
            ub      = [inf,inf,inf];

            weighted_deviations = @(p) (fitfun(p,histBinCenters(fitMask))-histValues(fitMask));
            optio = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','MaxFunctionEvaluations',100);
            try 
                [pFit,~,residFit,~,~,~,J]=lsqnonlin(weighted_deviations,pGuess,lb,ub,optio);
            catch
                display('fitting failed');
            end

            if plotDebug
                figure(889)
                subplot(1,2,dataOrControl)
                hold on
                bar(histBinCenters,histValues/max(histValues))
                plot([leftEdge,leftEdge],[0,1],'k--','LineWidth',2)
                plot([rightEdge,rightEdge],[0,1],'k--','LineWidth',2)
                plot(histBinCenters,fitfun(pFit,histBinCenters)/max(fitfun(pFit,histBinCenters)),'LineWidth',2)
                plot([pFit(2)+pFit(3),pFit(2)+pFit(3)],[0,1],'k-.','LineWidth',2)
                plot([pFit(2)-pFit(3),pFit(2)-pFit(3)],[0,1],'k-.','LineWidth',2)
                hold off
                xlim(10*[-pFit(3),pFit(3)])
                xlabel('Signal amplitude (V)')
                ylabel('rel. occourance (arb. units)')
                box on
                set(gca, 'FontName', 'Arial')
            set(gca,'FontSize', 12);
            end

            % find peaks
            dataToAnalyzeCentered = dataToAnalyze-pFit(2);
            
            snThreshold = pFit(3);
            
            [negSpikeHeigths,negSpikePos,negSpikeWidth,p]   = findpeaks(-dataToAnalyzeCentered,timesInSecondsToAnalyze,'MinPeakDistance',1/maxFireRate,'MinPeakProminence',SNThresholdNegative*pFit(3));
            [posSpikeHeigths,posSpikePos,posSpikeWidth,p2]  = findpeaks(dataToAnalyzeCentered,timesInSecondsToAnalyze,'MinPeakDistance',1/maxFireRate,'MinPeakProminence',SNThresholdPositive*pFit(3));

            if plotDebug
                figure(831)
                subplot(2,1,dataOrControl)
                hold on
                plot(timesInSecondsToAnalyze/60,dataToAnalyzeCentered,'-','LineWidth',1.5)
                errorbar(negSpikePos/60,-negSpikeHeigths,zeros(size(negSpikeHeigths)),zeros(size(negSpikeHeigths)),...
                         (negSpikeWidth/2)/60,(negSpikeWidth/2)/60,'.','MarkerSize',16,'Capsize',0,'LineWidth',2)
                plot(posSpikePos/60,posSpikeHeigths,'.','MarkerSize',16)
                hold off
                ylim([-20*pFit(3),10*pFit(3)])
                xlim([timesInSecondsToAnalyze(1)/60,timesInSecondsToAnalyze(1)/60+2])
                title('spike detection without postselection')
                legend('SG filtered data','detected spikes')
                box on
                xlabel('time (min)')
                ylabel('mV')
                set(gca, 'FontName', 'Arial')
                set(gca,'FontSize', 12);

            end

            % post select peaks 
            fullSpikeAmps = [];
            for jdx = 1:length(negSpikePos)
                nextPosSpike        = find(posSpikePos>negSpikePos(jdx),1,'first');
                if isempty(nextPosSpike)
                    fullSpikeAmps(jdx)  = negSpikeHeigths(jdx);
                else
                    fullSpikeAmps(jdx)  = posSpikeHeigths(nextPosSpike)+negSpikeHeigths(jdx);
                end

            end

            % cap fullSpikeAmps at 0.5V
            tooLageMask = fullSpikeAmps>0.4;
            fullSpikeAmps(tooLageMask)=0.4;

            if dataOrControl == 1
                % Data postselection: try to find peaks higher than norm:
                figure(771),clf;
                dataHistogram       = histogram(fullSpikeAmps,75); % ,'BinWidth',5e-3
                histValues          = dataHistogram.Values;
                histBinCenters      = dataHistogram.BinEdges(1:end-1)+mean(diff(dataHistogram.BinEdges));
                

                % prefit central peak:
                %fit central noise with gaussian to autothreshold signal and noise
                %[maxHistValue,maxIdx] = max(histValues);
                [~,maxIdx] = max(histValues);
                maxHistValue = mean(maxk(histValues,5));
                
                dataMask = find(histBinCenters<histBinCenters(maxIdx) & histValues>=0.3*maxHistValue,1,'first');
                leftEdge = histBinCenters(dataMask);
                if isempty(leftEdge)
                   leftEdge = 1; 
                end
                dataMask = find(histBinCenters>histBinCenters(maxIdx) & histValues>=0.3*maxHistValue,1,'last');
                rightEdge= histBinCenters(dataMask);

                fitMask         = histBinCenters>=leftEdge & histBinCenters<=rightEdge;
                fitfun          = @(p,x) p(1).*exp(-1/2*((x-p(2))/p(3)).^2);

                [maxValue, maxPos] = max(histValues);

                lb      = [0,-inf,0];  
                pGuess  = [maxValue,histBinCenters(maxPos),range([leftEdge,rightEdge]/2)];  % Inital (guess) parameters
                ub      = [inf,inf,inf];

                weighted_deviations = @(p) (fitfun(p,histBinCenters(fitMask))-histValues(fitMask));
                optio = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','MaxFunctionEvaluations',100);
                [pFit,~,residFit,~,~,~,J]=lsqnonlin(weighted_deviations,pGuess,lb,ub,optio);

                % fit whole distribution with double gaussian
                partfun         = @(p,x) 2*p(3)./(1+exp(p(4).*(x-p(2))));
                asymGauss       = @(p,x) p(1)./partfun(p,x).*exp(-4*log(2)*(2*((x-p(2))./partfun(p,x)).^2));
                parameterNames  = {'amp','peakPos','sigma','asym'};
                fitfun          = @(p,x) p(1).*exp(-1/2*((x-p(2))/p(3)).^2) + ...
                                         asymGauss([p(4),p(5),p(6),p(7)],x);

                lb      = [0,-inf,0,0,pFit(2)+normalEEGsigma*pFit(3),0,-inf];
                pGuess  = [pFit(1),pFit(2),pFit(3),pFit(1)/300,pFit(2)+5*pFit(3),6*pFit(3),-30];  % Inital (guess) parameters
                ub      = [inf,inf,inf,inf,inf,inf,0];

                weighted_deviations = @(p) (fitfun(p,histBinCenters(:))-histValues(:));
                optio = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','MaxFunctionEvaluations',100);
                [pFit,~,residFit,~,~,~,J]=lsqnonlin(weighted_deviations,pGuess,lb,ub,optio);

                if plotting
                    figure(690)
                    %subplot(1,2,idx)
                    hold on
                    bar(histBinCenters,histValues)
                    plot(histBinCenters,fitfun(pFit,histBinCenters),'LineWidth',2)
                    plot(histBinCenters,asymGauss(pFit(4:end),histBinCenters),'LineWidth',2)
                    plot([pFit(2)+normalEEGsigma*pFit(3),pFit(2)+normalEEGsigma*pFit(3)],[0,1.5*max(fitfun(pFit,histBinCenters))],'k-.','LineWidth',2)
                    plot([pFit(2)-normalEEGsigma*pFit(3),pFit(2)-normalEEGsigma*pFit(3)],[0,1.5*max(fitfun(pFit,histBinCenters))],'k-.','LineWidth',2)
                    hold off
                    box on
                    title(['file: ' files(fish).name])
                    xlabel('Signal amplitude (V)')
                    ylabel('occourance')
                    box on
                    set(gca, 'FontName', 'Arial')

                    if strcmp(exportFormat,'eps')
                        if ispc
                            print('-depsc', '-tiff', '-r300', '-painters', [DataFolder analysisFolder '\detectedSpikeHist_ID_' num2str(fish) '_' files(fish).name '.eps']);
                        elseif ismac
                            print('-depsc', '-tiff', '-r300', '-painters', [DataFolder analysisFolder '/detectedSpikeHist_ID_' num2str(fish) '_' files(fish).name '.eps']);
                        end
                    else
                        if ispc
                            print('-dpng', '-r300', '-painters', [DataFolder analysisFolder '\detectedSpikeHist_ID_' num2str(fish) '_' files(fish).name '.png']);
                        elseif ismac
                            print('-dpng', '-r300', '-painters', [DataFolder analysisFolder '/detectedSpikeHist_ID_' num2str(fish) '_' files(fish).name '.png']);
                        end
                    end
                end
            end
            
            highLowThreshold = pFit(2)+normalEEGsigma*pFit(3);
            
            gammaPeakMask = fullSpikeAmps>(pFit(2)+normalEEGsigma*pFit(3));
            betaPeakMask = fullSpikeAmps<(pFit(2)+normalEEGsigma*pFit(3));
            simplifiedSpikeTrace = zeros(size(fullSpikeAmps));
            simplifiedSpikeTrace(gammaPeakMask)=2;
            simplifiedSpikeTrace(betaPeakMask)=1;
            
            if plotting
                figure(832)
                if plotInMinutes
                    subplot(2,1,dataOrControl)
                    hold on
                    plot(timesInSecondsToAnalyze/60,dataToAnalyzeCentered,'-','LineWidth',0.75)
                    errorbar(negSpikePos(gammaPeakMask)/60,-negSpikeHeigths(gammaPeakMask),zeros(size(negSpikeHeigths(gammaPeakMask))),zeros(size(negSpikeHeigths(gammaPeakMask))),...
                             (negSpikeWidth(gammaPeakMask)/2)/60,(negSpikeWidth(gammaPeakMask)/2)/60,'.','MarkerSize',16,'Capsize',0,'LineWidth',2)
                    errorbar(negSpikePos(betaPeakMask)/60,-negSpikeHeigths(betaPeakMask),zeros(size(negSpikeHeigths(betaPeakMask))),zeros(size(negSpikeHeigths(betaPeakMask))),...
                             (negSpikeWidth(betaPeakMask)/2)/60,(negSpikeWidth(betaPeakMask)/2)/60,'.','MarkerSize',8,'Capsize',0,'LineWidth',2)

                    hold off
                    ylim([max(-30*pFit(3),1.1*min(dataToAnalyzeCentered)),min(10*pFit(3),1.1*max(dataToAnalyzeCentered))])
                    xlim([timesInSecondsToAnalyze(1)/60+5,timesInSecondsToAnalyze(1)/60+10])
                    title('post selected spikes')
                    legend('SG filtered data','\gamma spikes','\beta spikes')
                    box on
                    xlabel('time (min)')
                    ylabel('mV')
                    set(gca, 'FontName', 'Arial')
                    set(gca,'FontSize', 12);
                else
                    subplot(2,1,dataOrControl)
                    hold on
                    plot(timesInSecondsToAnalyze,dataToAnalyzeCentered,'-','LineWidth',0.75)
                    errorbar(negSpikePos(gammaPeakMask),-negSpikeHeigths(gammaPeakMask),zeros(size(negSpikeHeigths(gammaPeakMask))),zeros(size(negSpikeHeigths(gammaPeakMask))),...
                             (negSpikeWidth(gammaPeakMask)/2),(negSpikeWidth(gammaPeakMask)/2),'.','MarkerSize',16,'Capsize',0,'LineWidth',2)
                    errorbar(negSpikePos(betaPeakMask),-negSpikeHeigths(betaPeakMask),zeros(size(negSpikeHeigths(betaPeakMask))),zeros(size(negSpikeHeigths(betaPeakMask))),...
                             (negSpikeWidth(betaPeakMask)/2),(negSpikeWidth(betaPeakMask)/2),'.','MarkerSize',8,'Capsize',0,'LineWidth',2)

                    hold off
                    ylim([max(-30*pFit(3),1.1*min(dataToAnalyzeCentered)),min(10*pFit(3),1.1*max(dataToAnalyzeCentered))])
                    xlim([timesInSecondsToAnalyze(1)+5*60,timesInSecondsToAnalyze(1)+10*60])
                    title('post selected spikes')
                    legend('SG filtered data','\gamma spikes','\beta spikes')
                    box on
                    xlabel('time (min)')
                    ylabel('mV')
                    set(gca, 'FontName', 'Arial')
                    set(gca,'FontSize', 12);
                end

                if strcmp(exportFormat,'eps')
                    if ispc
                        print('-depsc', '-tiff', '-r300', '-painters', [DataFolder analysisFolder '\sampleSipkeDetection_ID_' num2str(fish) '_' files(fish).name '.eps']);
                    elseif ismac
                        print('-depsc', '-tiff', '-r300', '-painters', [DataFolder analysisFolder '/sampleSipkeDetection_ID_' num2str(fish) '_' files(fish).name '.eps']);
                    end
                else
                    if ispc
                        print('-dpng', '-r300', '-painters', [DataFolder analysisFolder '\sampleSipkeDetection_ID_' num2str(fish) '_' files(fish).name '.png']);
                    elseif ismac
                        print('-dpng', '-r300', '-painters', [DataFolder analysisFolder '/sampleSipkeDetection_ID_' num2str(fish) '_' files(fish).name '.png']);
                    end
                end


            end
            
            % FFT prototyping
            if plotDebug
                figure(819)
                subplot(3,2,dataOrControl)
                
                    histogram(diff(negSpikePos(~tooLageMask)),'BinWidth',0.1)
                    hold on
                    histogram(diff(negSpikePos(betaPeakMask)),'BinWidth',0.1)
                    plot([burstThresholdBeta,burstThresholdBeta],ylim,'k--')
                    hold off
                    legend('all spikes','beta spikes')
                    xlim([0,5])
                    box on
                    xlabel('\Deltat (s)')
                    ylabel('occurrence')
                    set(gca, 'FontName', 'Arial')
                    set(gca,'FontSize', 12);
                subplot(3,2,dataOrControl+2)
                    histogram(diff(negSpikePos(gammaPeakMask)),'BinWidth',0.1)
                    hold on
                    plot([burstThresholdBeta,burstThresholdBeta],ylim,'k--')
                    hold off
                    legend('all spikes','beta spikes')
                    xlim([0,5])
                    box on
                    xlabel('\Deltat (s)')
                    ylabel('occurrence')
                    set(gca, 'FontName', 'Arial')
                    set(gca,'FontSize', 12);
                subplot(3,2,dataOrControl+4)     
                    histogram(1./diff(negSpikePos(~tooLageMask)),75)
                    hold on
                    histogram(1./diff(negSpikePos(gammaPeakMask)),75)
                    histogram(1./diff(negSpikePos(betaPeakMask)),75)
                    plot([1/burstThresholdBeta,1/burstThresholdBeta],ylim,'k--')
                    hold off
                    legend('gamma spikes')
                    xlim([0,10])
                    box on
                    xlabel('(\Deltat)^{-1} (Hz)')
                    ylabel('occurrence')
                    set(gca, 'FontName', 'Arial')
                    set(gca,'FontSize', 12);
            end
           

            % 2020-03-28 beta response rep rate analysis
            figure(77112);
            
            numDomRepRates = 0;
            if dataOrControl == 1
                repRateHisto    = histogram(1./diff(negSpikePos(betaPeakMask)),50);
                repRateFreq     = repRateHisto.BinEdges(1:end-1)+diff(repRateHisto.BinEdges);
                repRateAmp      = repRateHisto.Values;
                
                repRateAmpFiltered = sgolayfilt(repRateAmp,3,9);
                
                plot(repRateFreq,repRateAmp)
                [peakRateAmp,peakRateFreq,peakRateBandwidth,~]   = findpeaks(repRateAmp,repRateFreq,'MinPeakProminence',10,'MinPeakHeight',mean(max(repRateAmp,5))/10);
                hold on
                plot(repRateFreq,repRateAmpFiltered)
                plot([median(1./diff(negSpikePos(betaPeakMask))),median(1./diff(negSpikePos(betaPeakMask)))],ylim,'k--')
                errorbar(peakRateFreq,peakRateAmp,zeros(size(peakRateAmp)),zeros(size(peakRateAmp)),peakRateBandwidth,peakRateBandwidth,'.','MarkerSize',16,'Capsize',0,'LineWidth',2)
                hold off
                
                xlabel('frequency (Hz)')
                ylabel('occurrence')
                set(gca, 'FontName', 'Arial')
                set(gca,'FontSize', 12);
                
                numDomRepRates = length(peakRateFreq);
                
            end
            
            
            
            
            % 2021-04-03 new burst analysis
            
            gamma_burstSpikeFreq            = [];
            gamma_burstDuration             = [];
            gamma_burstSpikes               = [];
            gamma_numTotalBursts_temp       = 0;
            gamma_IsolatedSpikesDuration_temp    = 0;
            gamma_IsolatedSpikes_temp            = 0;
            
            beta_burstSpikeFreq             = [];
            beta_burstDuration              = [];
            beta_burstSpikes                = [];
            beta_numTotalBursts_temp        = 0;
            beta_IsolatedSpikesDuration_temp     = 0;
            beta_IsolatedSpikes_temp             = 0;
            
            currentBurstNum     = 0;
            currentBurstPos     = [];
            lastPartOfBurstFlag = 0;
            
            for jdx = 2:length(negSpikePos)
                
                diffTime = negSpikePos(jdx)-negSpikePos(jdx-1);
                
                if ~lastPartOfBurstFlag
                    if currentBurstNum >1
                        if currentBurstType==2
                            gamma_burstSpikeFreq(end+1)   = mean(diff(currentBurstPos));
                            gamma_burstDuration(end+1)    = currentBurstPos(end)-currentBurstPos(1)+min(burstThresholdGamma,diffTime);
                            gamma_burstSpikes(end+1)      = currentBurstNum;
                            gamma_numTotalBursts_temp     = gamma_numTotalBursts_temp+1;
                        elseif currentBurstType==1 
                            beta_burstSpikeFreq(end+1)   = mean(diff(currentBurstPos));
                            beta_burstDuration(end+1)    = currentBurstPos(end)-currentBurstPos(1)+min(burstThresholdBeta,diffTime);
                            beta_burstSpikes(end+1)      = currentBurstNum;
                            beta_numTotalBursts_temp     = beta_numTotalBursts_temp+1;
                        end
                    elseif currentBurstNum == 1
                        if currentBurstType==2
                            gamma_IsolatedSpikesDuration_temp = gamma_IsolatedSpikesDuration_temp+min(burstThresholdGamma,diffTime);
                            gamma_IsolatedSpikes_temp = gamma_IsolatedSpikes_temp+1;
                        elseif currentBurstType==1
                            beta_IsolatedSpikesDuration_temp = beta_IsolatedSpikesDuration_temp+min(burstThresholdBeta,diffTime);
                            beta_IsolatedSpikes_temp = beta_IsolatedSpikes_temp+1;
                        end
                    end
                    currentBurstType = simplifiedSpikeTrace(jdx-1);
                    currentBurstNum = 1;
                    currentBurstPos = negSpikePos(jdx-1);
                    lastPartOfBurstFlag = 1;
                end
                    
                if(simplifiedSpikeTrace(jdx)==2 && currentBurstType==2)
                    if(diffTime<burstThresholdGamma)
                        currentBurstPos(end+1)  = negSpikePos(jdx);
                        currentBurstNum         = currentBurstNum+1;
                        lastPartOfBurstFlag     = 1;
                    else
                        lastPartOfBurstFlag     = 0;
                    end
                    
                elseif(simplifiedSpikeTrace(jdx)==1 && currentBurstType==1)
                    if(diffTime<burstThresholdGamma)
                        currentBurstPos(end+1)  = negSpikePos(jdx);
                        currentBurstNum         = currentBurstNum+1;
                        lastPartOfBurstFlag     = 1;
                    else
                        lastPartOfBurstFlag     = 0;
                    end
                    
                elseif(simplifiedSpikeTrace(jdx)==1 && currentBurstType==2)
                    lastPartOfBurstFlag     = 0;
                elseif(simplifiedSpikeTrace(jdx)==2 && currentBurstType==1)
                    lastPartOfBurstFlag     = 0;
                end
                
            end
           
            
            gamma_totalBurstDuration(fish,dataOrControl)  = sum(gamma_burstDuration);
            gamma_avBurstSpikeFreq(fish,dataOrControl)    = mean(gamma_burstSpikeFreq);
            gamma_avBurstDuration(fish,dataOrControl)     = mean(gamma_burstDuration);
            gamma_avBurstSpikes(fish,dataOrControl)       = mean(gamma_burstSpikes);
            gamma_medBurstSpikeFreq(fish,dataOrControl)    = median(gamma_burstSpikeFreq);
            gamma_medBurstDuration(fish,dataOrControl)     = median(gamma_burstDuration);
            gamma_medBurstSpikes(fish,dataOrControl)       = median(gamma_burstSpikes);
            gamma_numTotalBursts(fish,dataOrControl)      = gamma_numTotalBursts_temp;
            gamma_burstFreq(fish,dataOrControl)           = gamma_numTotalBursts_temp./(timesInSecondsToAnalyze(end)-timesInSecondsToAnalyze(1));
            
            gamma_IsolatedSpikesDuration(fish,dataOrControl) = gamma_IsolatedSpikesDuration_temp;
            gamma_IsolatedSpikes(fish,dataOrControl)         = gamma_IsolatedSpikes_temp;
            
            gamma_spikeRateInHz(fish,dataOrControl)     = length(negSpikeHeigths(gammaPeakMask))/(timesInSecondsToAnalyze(end)-timesInSecondsToAnalyze(1));
            gamma_spikeStrength(fish,dataOrControl)     = mean(negSpikeHeigths(gammaPeakMask)'.*negSpikeWidth(gammaPeakMask));
            gamma_spikeDuration(fish,dataOrControl)     = mean(negSpikeWidth(gammaPeakMask));
            gamma_spikeAmplitude(fish,dataOrControl)    = mean(negSpikeHeigths(gammaPeakMask));
            
            gamma_spikeRateInHz_med(fish,dataOrControl)     = median(1./diff(negSpikePos(gammaPeakMask)));
            gamma_spikeStrength_med(fish,dataOrControl)     = median(negSpikeHeigths(gammaPeakMask)'.*negSpikeWidth(gammaPeakMask));
            gamma_spikeDuration_med(fish,dataOrControl)     = median(negSpikeWidth(gammaPeakMask));
            gamma_spikeAmplitude_med(fish,dataOrControl)    = median(negSpikeHeigths(gammaPeakMask));
            
            beta_totalBurstDuration(fish,dataOrControl)  = sum(beta_burstDuration);
            beta_avBurstSpikeFreq(fish,dataOrControl)    = mean(beta_burstSpikeFreq);
            beta_avBurstDuration(fish,dataOrControl)     = mean(beta_burstDuration);
            beta_avBurstSpikes(fish,dataOrControl)       = mean(beta_burstSpikes);
            beta_medBurstSpikeFreq(fish,dataOrControl)    = median(beta_burstSpikeFreq);
            beta_medBurstDuration(fish,dataOrControl)     = median(beta_burstDuration);
            beta_medBurstSpikes(fish,dataOrControl)       = median(beta_burstSpikes);
            beta_numTotalBursts(fish,dataOrControl)      = beta_numTotalBursts_temp;
            beta_burstFreq(fish,dataOrControl)           = beta_numTotalBursts_temp./(timesInSecondsToAnalyze(end)-timesInSecondsToAnalyze(1));
            
            beta_IsolatedSpikesDuration(fish,dataOrControl) = beta_IsolatedSpikesDuration_temp;
            beta_IsolatedSpikes(fish,dataOrControl)         = beta_IsolatedSpikes_temp;
            
            beta_spikeRateInHz(fish,dataOrControl)       = length(negSpikeHeigths(betaPeakMask))/(timesInSecondsToAnalyze(end)-timesInSecondsToAnalyze(1));
            beta_spikeStrength(fish,dataOrControl)       = mean(negSpikeHeigths(betaPeakMask)'.*negSpikeWidth(betaPeakMask));
            beta_spikeDuration(fish,dataOrControl)       = mean(negSpikeWidth(betaPeakMask));
            beta_spikeAmplitude(fish,dataOrControl)      = mean(negSpikeHeigths(betaPeakMask));
            
            beta_spikeRateInHz_med(fish,dataOrControl)       = median(1./diff(negSpikePos(betaPeakMask)));
            beta_spikeStrength_med(fish,dataOrControl)       = median(negSpikeHeigths(betaPeakMask)'.*negSpikeWidth(betaPeakMask));
            beta_spikeDuration_med(fish,dataOrControl)       = median(negSpikeWidth(betaPeakMask));
            beta_spikeAmplitude_med(fish,dataOrControl)      = median(negSpikeHeigths(betaPeakMask));
            
            total_time                                   = (timesInSecondsToAnalyze(end)-timesInSecondsToAnalyze(1));
            ratio_beta(fish,dataOrControl)               = beta_totalBurstDuration(fish,dataOrControl)./total_time+beta_IsolatedSpikesDuration(fish,dataOrControl)./total_time;
            ratio_gamma(fish,dataOrControl)              = gamma_totalBurstDuration(fish,dataOrControl)./total_time+gamma_IsolatedSpikesDuration(fish,dataOrControl)./total_time;
            ratio_alpha(fish,dataOrControl)              = 1-ratio_beta(fish,dataOrControl)-ratio_gamma(fish,dataOrControl) ;
            
            % old burst analysis (will be removed in future)
            
            % 2021-04-11 gamma bursts on top of beta analysis
            selectedSpikePos    = negSpikePos(gammaPeakMask);
            currentBurstNum     = 0;
            currentBurstPos     = [];
            gamma_ignore_beta_burstSpikeFreq      = [];
            gamma_ignore_beta_burstDuration       = [];
            gamma_ignore_beta_burstSpikes         = [];
            gamma_ignore_beta_numTotalBursts_temp = 0;
            lastPartOfBurstFlag = 0;

            for jdx = 2:length(selectedSpikePos)
                diffTime = selectedSpikePos(jdx)-selectedSpikePos(jdx-1);
                if(diffTime<burstThresholdGamma)
                    currentBurstPos(end+1)  = negSpikePos(jdx);
                    currentBurstNum         = currentBurstNum+1;
                    lastPartOfBurstFlag     = 1;
                else
                    if lastPartOfBurstFlag && currentBurstNum >1
                        gamma_ignore_beta_burstSpikeFreq(end+1)   = mean(diff(currentBurstPos));
                        gamma_ignore_beta_burstDuration(end+1)    = currentBurstPos(end)-currentBurstPos(1)+burstThresholdGamma;
                        gamma_ignore_beta_burstSpikes(end+1)      = currentBurstNum;
                        gamma_ignore_beta_numTotalBursts_temp     = gamma_ignore_beta_numTotalBursts_temp+1;
                    end
                    currentBurstNum = 0;
                    currentBurstPos = [];
                    lastPartOfBurstFlag = 0;
                end
            end

            gamma_ignore_beta_totalBurstDuration(fish,dataOrControl)  = sum(gamma_ignore_beta_burstDuration);
            gamma_ignore_beta_avBurstSpikeFreq(fish,dataOrControl)    = mean(gamma_ignore_beta_burstSpikeFreq);
            gamma_ignore_beta_avBurstDuration(fish,dataOrControl)     = mean(gamma_ignore_beta_burstDuration);
            gamma_ignore_beta_avBurstSpikes(fish,dataOrControl)       = mean(gamma_ignore_beta_burstSpikes);
            gamma_ignore_beta_numTotalBursts(fish,dataOrControl)      = gamma_ignore_beta_numTotalBursts_temp;
            gamma_ignore_beta_burstFreq(fish,dataOrControl)           = gamma_ignore_beta_numTotalBursts_temp./(timesInSecondsToAnalyze(end)-timesInSecondsToAnalyze(1));
            
           
            %export data
            
            
            
            if dataOrControl == 1
                myFishData.fish             = fish;
                myFishData.filename         = {files(fish).name};
                myFishData.block = ['first ' num2str(startTimeBlockInMinutes,'%03.f') ' min'];
                
                myFishData.ratio_alpha                = ratio_alpha(fish,dataOrControl);
                myFishData.ratio_beta                 = ratio_beta(fish,dataOrControl);
                myFishData.ratio_gamma                = ratio_gamma(fish,dataOrControl);

                myFishData.gamma_spikeRateInHz        = gamma_spikeRateInHz(fish,dataOrControl);
                myFishData.gamma_spikeStrength        = gamma_spikeStrength(fish,dataOrControl);
                myFishData.gamma_spikeDuration        = gamma_spikeDuration(fish,dataOrControl);
                myFishData.gamma_spikeAmplitude       = gamma_spikeAmplitude(fish,dataOrControl);
                
                myFishData.gamma_spikeRateInHz_med        = gamma_spikeRateInHz_med(fish,dataOrControl);
                myFishData.gamma_spikeStrength_med        = gamma_spikeStrength_med(fish,dataOrControl);
                myFishData.gamma_spikeDuration_med        = gamma_spikeDuration_med(fish,dataOrControl);
                myFishData.gamma_spikeAmplitude_med       = gamma_spikeAmplitude_med(fish,dataOrControl);
                
                myFishData.gamma_avBurstSpikeFreq     = gamma_avBurstSpikeFreq(fish,dataOrControl);
                myFishData.gamma_avBurstDuration      = gamma_avBurstDuration(fish,dataOrControl);
                myFishData.gamma_avBurstSpikes        = gamma_avBurstSpikes(fish,dataOrControl);
                myFishData.gamma_medBurstSpikeFreq     = gamma_medBurstSpikeFreq(fish,dataOrControl);
                myFishData.gamma_medBurstDuration      = gamma_medBurstDuration(fish,dataOrControl);
                myFishData.gamma_medBurstSpikes        = gamma_medBurstSpikes(fish,dataOrControl);
                myFishData.gamma_totalBurstDuration   = gamma_totalBurstDuration(fish,dataOrControl);
                myFishData.gamma_numTotalBursts       = gamma_numTotalBursts(fish,dataOrControl);
                myFishData.gamma_burstFreq            = gamma_burstFreq(fish,dataOrControl);
                
                myFishData.gamma_IsolatedSpikesDuration = gamma_IsolatedSpikesDuration(fish,dataOrControl);
                myFishData.gamma_IsolatedSpikes         = gamma_IsolatedSpikes(fish,dataOrControl);
            

                myFishData.beta_spikeRateInHz         = beta_spikeRateInHz(fish,dataOrControl);
                myFishData.beta_spikeStrength         = beta_spikeStrength(fish,dataOrControl);
                myFishData.beta_spikeDuration         = beta_spikeDuration(fish,dataOrControl);
                myFishData.beta_spikeAmplitude        = beta_spikeAmplitude(fish,dataOrControl);
                
                myFishData.beta_spikeRateInHz_med         = beta_spikeRateInHz_med(fish,dataOrControl);
                myFishData.beta_spikeStrength_med         = beta_spikeStrength_med(fish,dataOrControl);
                myFishData.beta_spikeDuration_med         = beta_spikeDuration_med(fish,dataOrControl);
                myFishData.beta_spikeAmplitude_med        = beta_spikeAmplitude_med(fish,dataOrControl);

                myFishData.beta_avBurstSpikeFreq      = beta_avBurstSpikeFreq(fish,dataOrControl);
                myFishData.beta_avBurstDuration       = beta_avBurstDuration(fish,dataOrControl);
                myFishData.beta_avBurstSpikes         = beta_avBurstSpikes(fish,dataOrControl);
                myFishData.beta_medBurstSpikeFreq      = beta_medBurstSpikeFreq(fish,dataOrControl);
                myFishData.beta_medBurstDuration       = beta_medBurstDuration(fish,dataOrControl);
                myFishData.beta_medBurstSpikes         = beta_medBurstSpikes(fish,dataOrControl);
                myFishData.beta_totalBurstDuration    = beta_totalBurstDuration(fish,dataOrControl);
                myFishData.beta_numTotalBursts        = beta_numTotalBursts(fish,dataOrControl);
                myFishData.beta_burstFreq             = beta_burstFreq(fish,dataOrControl);
                
                myFishData.beta_IsolatedSpikesDuration = beta_IsolatedSpikesDuration(fish,dataOrControl);
                myFishData.beta_IsolatedSpikes         = beta_IsolatedSpikes(fish,dataOrControl);
                
                
                myFishData.gamma_ignore_beta_totalBurstDuration     = gamma_ignore_beta_totalBurstDuration(fish,dataOrControl); 
                myFishData.gamma_ignore_beta_avBurstSpikeFreq     = gamma_ignore_beta_avBurstSpikeFreq(fish,dataOrControl);   
                myFishData.gamma_ignore_beta_avBurstDuration     = gamma_ignore_beta_avBurstDuration(fish,dataOrControl);    
                myFishData.gamma_ignore_beta_avBurstSpikes     = gamma_ignore_beta_avBurstSpikes(fish,dataOrControl);
                myFishData.gamma_ignore_beta_numTotalBursts     = gamma_ignore_beta_numTotalBursts(fish,dataOrControl);
                myFishData.gamma_ignore_beta_burstFreq     = gamma_ignore_beta_burstFreq(fish,dataOrControl);



                myFishData.beta_numDomRepRates        = numDomRepRates;
                myFishData.snThreshold                = snThreshold;
                myFishData.highLowThreshold           = highLowThreshold;

                
                
            else
                myFishDataEnd.fish             = fish;
                myFishDataEnd.filename         = {files(fish).name};
                myFishDataEnd.block = ['last  ' num2str(endTimeBlockInMinutes,'%03.f') ' min'];
                
                myFishDataEnd.ratio_alpha                = ratio_alpha(fish,dataOrControl);
                myFishDataEnd.ratio_beta                 = ratio_beta(fish,dataOrControl);
                myFishDataEnd.ratio_gamma                = ratio_gamma(fish,dataOrControl);

                myFishDataEnd.gamma_spikeRateInHz        = gamma_spikeRateInHz(fish,dataOrControl);
                myFishDataEnd.gamma_spikeStrength        = gamma_spikeStrength(fish,dataOrControl);
                myFishDataEnd.gamma_spikeDuration        = gamma_spikeDuration(fish,dataOrControl);
                myFishDataEnd.gamma_spikeAmplitude       = gamma_spikeAmplitude(fish,dataOrControl);
                
                myFishDataEnd.gamma_avBurstSpikeFreq     = gamma_avBurstSpikeFreq(fish,dataOrControl);
                myFishDataEnd.gamma_avBurstDuration      = gamma_avBurstDuration(fish,dataOrControl);
                myFishDataEnd.gamma_avBurstSpikes        = gamma_avBurstSpikes(fish,dataOrControl);
                myFishDataEnd.gamma_totalBurstDuration   = gamma_totalBurstDuration(fish,dataOrControl);
                myFishDataEnd.gamma_numTotalBursts       = gamma_numTotalBursts(fish,dataOrControl);
                myFishDataEnd.gamma_burstFreq            = gamma_burstFreq(fish,dataOrControl);
                
                myFishDataEnd.gamma_IsolatedSpikesDuration = gamma_IsolatedSpikesDuration(fish,dataOrControl);
                myFishDataEnd.gamma_IsolatedSpikes         = gamma_IsolatedSpikes(fish,dataOrControl);
            

                myFishDataEnd.beta_spikeRateInHz         = beta_spikeRateInHz(fish,dataOrControl);
                myFishDataEnd.beta_spikeStrength         = beta_spikeStrength(fish,dataOrControl);
                myFishDataEnd.beta_spikeDuration         = beta_spikeDuration(fish,dataOrControl);
                myFishDataEnd.beta_spikeAmplitude        = beta_spikeAmplitude(fish,dataOrControl);

                myFishDataEnd.beta_avBurstSpikeFreq      = beta_avBurstSpikeFreq(fish,dataOrControl);
                myFishDataEnd.beta_avBurstDuration       = beta_avBurstDuration(fish,dataOrControl);
                myFishDataEnd.beta_avBurstSpikes         = beta_avBurstSpikes(fish,dataOrControl);
                myFishDataEnd.beta_totalBurstDuration    = beta_totalBurstDuration(fish,dataOrControl);
                myFishDataEnd.beta_numTotalBursts        = beta_numTotalBursts(fish,dataOrControl);
                myFishDataEnd.beta_burstFreq             = beta_burstFreq(fish,dataOrControl);
                
                myFishDataEnd.beta_IsolatedSpikesDuration = beta_IsolatedSpikesDuration(fish,dataOrControl);
                myFishDataEnd.beta_IsolatedSpikes         = beta_IsolatedSpikes(fish,dataOrControl);


                myFishDataEnd.beta_numDomRepRates        = numDomRepRates;
                myFishDataEnd.snThreshold                = snThreshold;
                myFishDataEnd.highLowThreshold           = highLowThreshold;

                myFishDataEnd.gamma_ignore_beta_totalBurstDuration     = gamma_ignore_beta_totalBurstDuration(fish,dataOrControl); 
                myFishDataEnd.gamma_ignore_beta_avBurstSpikeFreq     = gamma_ignore_beta_avBurstSpikeFreq(fish,dataOrControl);   
                myFishDataEnd.gamma_ignore_beta_avBurstDuration     = gamma_ignore_beta_avBurstDuration(fish,dataOrControl);    
                myFishDataEnd.gamma_ignore_beta_avBurstSpikes     = gamma_ignore_beta_avBurstSpikes(fish,dataOrControl);
                myFishDataEnd.gamma_ignore_beta_numTotalBursts     = gamma_ignore_beta_numTotalBursts(fish,dataOrControl);
                myFishDataEnd.gamma_ignore_beta_burstFreq     = gamma_ignore_beta_burstFreq(fish,dataOrControl);
                
            end
            
            
            
            myAnalysisSettings = struct;
            myAnalysisSettings.ExcludeInitialMinutes   = ExcludeInitialMinutes;
            myAnalysisSettings.ExcludeLastMinutes      = ExcludeLastMinutes;
            myAnalysisSettings.startTimeBlockInMinutes = startTimeBlockInMinutes;
            myAnalysisSettings.endTimeBlockInMinutes   = endTimeBlockInMinutes;
            myAnalysisSettings.SNThresholdNegative     = SNThresholdNegative;
            myAnalysisSettings.SNThresholdPositive     = SNThresholdPositive;
            myAnalysisSettings.maxFireRate             = maxFireRate;
            myAnalysisSettings.downsamplingFactor      = downsamplingFactor;
            myAnalysisSettings.alpha                   = alpha;
            myAnalysisSettings.SGorder                 = SGorder;
            myAnalysisSettings.SGlength                = SGlength;
            myAnalysisSettings.normalEEGsigma          = normalEEGsigma;
            myAnalysisSettings.plotting                = plotting;
            myAnalysisSettings.plotDebug               = plotDebug;
            myAnalysisSettings.plotInMinutes           = plotInMinutes;
            myAnalysisSettings.exportFormat            = exportFormat;
            
            
            
            drawnow
            

        end
        
        exportTable     = [exportTable; myFishData];
        exportTableEnd  = [exportTableEnd; myFishDataEnd];
        
        if pauseEachFile
            figure(111),clf;
            waitforbuttonpress
        end
    end

    if ispc
        writetable(exportTable, [DataFolder analysisFolder '\' excelTableName '.csv']);
        writetable(exportTableEnd, [DataFolder analysisFolder '\' excelTableName '_end.csv']);
        save([DataFolder analysisFolder '\myFishData' ],'exportTable')
        save([DataFolder analysisFolder '\myFishData_end' ],'exportTableEnd')
        save([DataFolder analysisFolder '\myAnalysisSettings' ],'myAnalysisSettings')
    elseif ismac
        writetable(exportTable, [DataFolder analysisFolder '/' excelTableName '.csv']);
        writetable(exportTableEnd, [DataFolder analysisFolder '/' excelTableName '_end.csv']);
        save([DataFolder analysisFolder '/myFishData' ],'exportTable')
        save([DataFolder analysisFolder '/myFishData_end' ],'exportTableEnd')
        save([DataFolder analysisFolder '/myAnalysisSettings' ],'myAnalysisSettings')
    end


    % TODO: 
    % format output
    % some buttons 

end
