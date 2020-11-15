function EEG_SpikeFinder(DataFolder,varargin)
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

    % plotting settings
    p.addParameter('plotting',                  true);  % 
    p.addParameter('plotDebug',                 true);  % 
    p.addParameter('plotInMinutes',             true);  %

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
    
    plotting                = p.Results.plotting;
    plotDebug               = p.Results.plotDebug;
    plotInMinutes           = p.Results.plotInMinutes;
    
    exportFormat            = p.Results.exportFormat;
   
    
    
    % create empty containers
    spikeRateInHz           = [];
    spikeStrength           = [];
    spikeDuration           = [];
    spikeAmplitude          = [];
    exportTable             = table;
    excelTableName          = 'EEG_SpikeFinder_results';

    %make new analysis folder
    formatOut       = 'yyyy-mm-dd_hh-MM';
    timeString      = datestr(now,formatOut);
    analysisFolder  = [timeString '_analysis'];

    mkdir([DataFolder analysisFolder])

    files = dir([DataFolder '*.abf']);
    for fish=1:length(files)
        myFishData = table;
        % load data
        [data, si, header]       = abfload([DataFolder files(fish).name]);
        % get sampling frequency
        samplingFrequency   = 1/header.si*1e6; % in Hz 

        % downsample data
        dataDownSampled = imresize(data,1/downsamplingFactor,'box');
        timesInSeconds = 1/samplingFrequency*(0:downsamplingFactor:length(data)-1);

        % run SG filter over data
        dataSGfiltered  = sgolayfilt(dataDownSampled,SGorder,SGlength);

        if plotDebug
            figure(881),clf;

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

            fitMask         = histBinCenters>=leftEdge & histBinCenters<=rightEdge;
            fitfun          = @(p,x) p(1).*exp(-1/2*((x-p(2))/p(3)).^2);

            [maxValue, maxPos] = max(histValues);

            lb      = [0,-inf,0];  
            pGuess  = [maxValue,histBinCenters(maxPos),range([leftEdge,rightEdge]/2)];  % Inital (guess) parameters
            ub      = [inf,inf,inf];

            weighted_deviations = @(p) (fitfun(p,histBinCenters(fitMask))-histValues(fitMask));
            optio = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','MaxFunctionEvaluations',100);
            [pFit,~,residFit,~,~,~,J]=lsqnonlin(weighted_deviations,pGuess,lb,ub,optio);

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
                [maxValue, maxPos]  = max(histValues);

                % prefit central peak:
                %fit central noise with gaussian to autothreshold signal and noise
                [maxHistValue,maxIdx] = max(histValues);
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

                [maxValue, maxPos] = max(histValues);

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
            peakMask = fullSpikeAmps>(pFit(2)+normalEEGsigma*pFit(3));

            if plotting
                figure(832)
                subplot(2,1,dataOrControl)
                hold on
                plot(timesInSecondsToAnalyze/60,dataToAnalyzeCentered,'-','LineWidth',0.75)
                errorbar(negSpikePos(peakMask)/60,-negSpikeHeigths(peakMask),zeros(size(negSpikeHeigths(peakMask))),zeros(size(negSpikeHeigths(peakMask))),...
                         (negSpikeWidth(peakMask)/2)/60,(negSpikeWidth(peakMask)/2)/60,'.','MarkerSize',16,'Capsize',0,'LineWidth',2)

                hold off
                ylim([max(-30*pFit(3),1.1*min(dataToAnalyzeCentered)),min(10*pFit(3),1.1*max(dataToAnalyzeCentered))])
                xlim([timesInSecondsToAnalyze(1)/60+5,timesInSecondsToAnalyze(1)/60+10])
                title('post selected spikes')
                legend('SG filtered data','detected spikes')
                box on
                xlabel('time (min)')
                ylabel('mV')
                set(gca, 'FontName', 'Arial')
                set(gca,'FontSize', 12);

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


            % save data
            spikeRateInHz(fish,dataOrControl)     = length(negSpikeHeigths(peakMask))/(timesInSecondsToAnalyze(end)-timesInSecondsToAnalyze(1));
            spikeStrength(fish,dataOrControl)     = mean(negSpikeHeigths(peakMask)'.*negSpikeWidth(peakMask));
            spikeDuration(fish,dataOrControl)     = mean(negSpikeWidth(peakMask));
            spikeAmplitude(fish,dataOrControl)    = mean(negSpikeHeigths(peakMask));


            % burst analysis 2020-10-21
            selectedSpikePos    = negSpikePos(peakMask);
            currentBurstNum     = 0;
            currentBurstPos     = [];
            burstSpikeFreq      = [];
            burstDuration       = [];
            burstSpikes         = [];
            burstThreshold      = 1;
            numTotalBursts_temp      = 0;
            lastPartOfBurstFlag = 0;

            for jdx = 2:length(selectedSpikePos)
                diffTime = selectedSpikePos(jdx)-selectedSpikePos(jdx-1);
                if(diffTime<burstThreshold)
                    currentBurstPos(end+1)  = negSpikePos(jdx);
                    currentBurstNum         = currentBurstNum+1;
                    lastPartOfBurstFlag     = 1;
                else
                    if lastPartOfBurstFlag && currentBurstNum >1
                        burstSpikeFreq(end+1)   = mean(diff(currentBurstPos));
                        burstDuration(end+1)    = currentBurstPos(end)-currentBurstPos(1);
                        burstSpikes(end+1)      = currentBurstNum;
                        numTotalBursts_temp          = numTotalBursts_temp+1;
                    end
                    currentBurstNum = 0;
                    currentBurstPos = [];
                    lastPartOfBurstFlag = 0;
                end
            end

            avBurstSpikeFreq(fish,dataOrControl)    = mean(burstSpikeFreq);
            avBurstDuration(fish,dataOrControl)     = mean(burstDuration);
            avBurstSpikes(fish,dataOrControl)       = mean(burstSpikes);
            numTotalBursts(fish,dataOrControl)      = numTotalBursts_temp;
            burstFreq(fish,dataOrControl)           = numTotalBursts_temp./(timesInSecondsToAnalyze(end)-timesInSecondsToAnalyze(1));

            %export data
            myFishData.fish             = fish;
            if dataOrControl == 1
                myFishData.block = ['first ' num2str(startTimeBlockInMinutes,'%03.f') ' min'];
            else
                myFishData.block = ['last  ' num2str(endTimeBlockInMinutes,'%03.f') ' min'];
            end
            myFishData.spikeRateInHz        = spikeRateInHz(fish,dataOrControl);
            myFishData.spikeStrength        = spikeStrength(fish,dataOrControl);
            myFishData.spikeDuration        = spikeDuration(fish,dataOrControl);
            myFishData.spikeAmplitude       = spikeAmplitude(fish,dataOrControl);

            myFishData.avBurstSpikeFreq     = avBurstSpikeFreq(fish,dataOrControl);
            myFishData.avBurstDuration      = avBurstDuration(fish,dataOrControl);
            myFishData.avBurstSpikes        = avBurstSpikes(fish,dataOrControl);
            myFishData.numTotalBursts       = numTotalBursts(fish,dataOrControl);
            myFishData.burstFreq            = burstFreq(fish,dataOrControl);

            myFishData.filename             = {files(fish).name};
            
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
            
            exportTable = [exportTable; myFishData];
            
            drawnow

        end

    end

    if ispc
        writetable(exportTable, [DataFolder analysisFolder '\' excelTableName '.csv']);
        save([DataFolder analysisFolder '\myFishData' ],'myFishData')
        save([DataFolder analysisFolder '\myAnalysisSettings' ],'myAnalysisSettings')
    elseif ismac
        writetable(exportTable, [DataFolder analysisFolder '/' excelTableName '.csv']);
        save([DataFolder analysisFolder '/myFishData' ],'myFishData')
        save([DataFolder analysisFolder '/myAnalysisSettings' ],'myAnalysisSettings')
    end


    fprintf('-------------------------------------\n');
    fprintf('mean spike rate (Hz): %.2f +/- %.2f \n',mean(spikeRateInHz(:,1),1),std(spikeRateInHz(:,1),1)./length(spikeRateInHz(:,1)));
    fprintf('-------------------------------------\n');

end
