classdef SignalAnalysis
    %Simplified SEE SignalanalysisClass
    
    properties
        %Filename can be the pwd address
        Filename
        %the signal for testing, vary with the resample,denoise
        Test_signal
        %test time
        Test_time
        %test frequency
        Test_freq
        %Background Subtraction
        Test_signal_offset
        %Backgorund by rloess method
        Test_signal_trendline
        %Reaction satates specify as 'Oxi' or 'Red'
        Reactionstates
        %the 2 times std of offset back
        BackGroundHeightThr
        %Flip find the signal interval list as index
        FlipFindInte
        %found peaks properties height location width(
        PeaksFound
        %Trainset Raw Tem
        TrainingSignalOri
        %mean value for the raw templates clustering
        RawTemFeaturesMu
        %varence value for the raw templates clustering
        RawTemFeaturesSigma
        %Raw templates clustering features (4 basic features +  height
        %width relative location)
        RawTemClusterFeatures
        %Feature normolized
        RawTemClusterFeaturesNorm
        %Numbers of raw unregulated templates
        NumTemplateRawUnRegu
        %Unregulared Raw Unregulated templates
        TemplateRawUnRegu
        %Counts for each Unregulated templates
        TemplateRawUnReguCounts
        %Regulated Raw templates
        TemplatateRawRegu
        %Regulated Raw templates Num
        TemplatateRawReguNum
        %Similarity cell
        SimilarityLag
        %Template Matched interval after filtering
        TemplatateMatchedInte
        %Template Matched interval Signal
        AMInteSig
        %After template matching traning set
        AMTraningSet
        %After template matching traning set Norm
        AMTraningSetNorm
        %After template matching traning set Sigma
        AMTemPSigma
        %After template matching traning set Mu
        AMTemPMu
        %After template matching ReNormed centroids
        AMTemCen
        %After template matching Regenerate templates
        AMTemplates
        %After matching the interval index
        AMTemidx
%         %Data comuted by conventional method by computing area of Guassian
%         %distribution
%         DataGaussain
%         %simulated gaussain siganl
%         SigGaussain
        
    end
    properties(Hidden,Access = private)
        gpuStates
        BinaryY
        Intecolors
        AMIntecolors
    end
    
    methods
        function thisSignalAnalysis = SignalAnalysis(Filename)
            %Input the type of files of the raw signal'*.txt';
            %thisSignalAnalysis = SignalAnalysis(Filename)
            % Input -- File name as a string, by defalut, consider there are other text string in the data file
            % Output -- Sig.test_signal Original signal vector in original unit
            %           Sig.Test_time Original time vector in original unit
            %           Sig.Test_freq test frequency depends on the test_time vector
            try
                thisSignalAnalysis.gpuStates = gpuDevice;
            catch
                thisSignalAnalysis.gpuStates = [];
            end
            thisSignalAnalysis.BinaryY = ["Yes","yes","YES","Y","y","OK","ok"];
            thisSignalAnalysis.Filename = Filename;
            Data = importdata(Filename).data;
            thisSignalAnalysis.Test_signal = Data(:,2);
            thisSignalAnalysis.Test_time = Data(:,1);
            thisSignalAnalysis.Test_freq  = 1/mean(diff(thisSignalAnalysis.Test_time));
            
        end
        function [theSignalAnalysis] = Preprocess(theSignalAnalysis,OriginalUnit,Timeshift,Resample,ResamplePlot)
            %[theSignalAnalysis] = Preprocess(theSignalAnalysis,OriginalUnit,Timeshift,Resample)
            % preprocess for the SEE signal, change the unit, shift the
            % time axis, and resample the data in case uneven sampling
%             Input -- OriginalUnit a string of orgiginal data's unit e.g. 'pA', it will convert different current into pA in this toolbox
%                      Timeshift wether or not performing the time shift, will shift the test_time starts from 0, input string as e.g. 'y'
%                      Resample wether or not to perform the resampling for uneven sampling, input string as e.g. 'y'
%                      ResamplePlot wether or not to show the plot with and without resampling, input string as e.g. 'y'
%             Output -- Sig.test_signal converted unit signal vector, depending on chosen option becomes resampled signal
%                       Sig.Test_time updated time vector in original unit, depending on chosen option becomes shifted time from 0
%                       Sig.Test_freq updatedtest frequency depends on the test_time vector
            if nargin == 1
                OriginalUnit = "A" ;
                Timeshift = "N";
                Resample = "N";
                ResamplePlot = "N";
            end
            switch OriginalUnit
                case "A"
                    Factor = 1e12;
                case {"mA","MA"}
                    Factor = 1e9;
                case {"uA","UA"}
                    Factor = 1e6;
                case {"nA","NA"}
                    Factor = 1e3;
                case {"FA","fA"}
                    Factor = 1e-3;
                case {"aA","AA"}
                    Factor = 1e-6;
                case {"PA","pA"}
                    Factor = 1;
                    
            end
            theSignalAnalysis.Test_signal = theSignalAnalysis.Test_signal * Factor;
            
            if ismember(Timeshift,theSignalAnalysis.BinaryY)
                theSignalAnalysis.Test_time = ...
                    [0:1/theSignalAnalysis.Test_freq:(length(theSignalAnalysis.Test_signal)-1)* 1/theSignalAnalysis.Test_freq]';
            end
            
            if ismember(Resample,theSignalAnalysis.BinaryY)
                diff_all = diff(theSignalAnalysis.Test_time);
                Ts = mode(diff_all);
                Fs = round(1/Ts);
                x = theSignalAnalysis.Test_signal;
                ts = theSignalAnalysis.Test_time - theSignalAnalysis.Test_time(1);
                y = resample(gather(x),gather(ts),gather(Fs));
                ts2 = (0:(length(y)-1))*Ts;
                if ismember(ResamplePlot,theSignalAnalysis.BinaryY)
                    figure;
                    plot(ts(10:end-15),x(10:end-15)); hold on;
                    plot(ts2(10:end-15),y(10:end-15));
                    title('Original Curve');
                    xlabel('Time [s]');
                    ylabel('Current [pA]');
                    %gpu arrary debug
                    minX = gather(min(ts(:,1)));
                    maxX = gather(max(ts(:,1)));
                    xlim([minX, maxX]);
                    legend('Raw data','Resampled');
                    legend('boxoff')
                end
                theSignalAnalysis.Test_signal = y(10:end-15);
                theSignalAnalysis.Test_time = ts2(10:end-15)';
                theSignalAnalysis.Test_freq = Fs;
            end
            
        end
        
        function [theSignalAnalysis] = Denoise(theSignalAnalysis,CutoffFreq,FilterOrder)
            %[theSignalAnalysis] = Denoise(theSignalAnalysis,CutoffFreq,FilterOrder)
            %Denoise signal bsed on FIR lsq low pass filter, changing filterorder to see the result
            %            Input --
            %            CutoffFreq the cut-off frequency to apply for FIR filter, the spike frequency can get from stft plot by leaving CutoffFreq empty, suggest running in .m since the in .mlx the plot won't show instantly e.g. 10
            %            FilterOrder the order for the FIR filter also known as term number,  higher the order, shaper the filter at pass band, but higher the delay, we din't perform the zero-phase filtering in this toolbox e.g. 10, can be changed in the command window it will ask wether or not to change, if enter no or nothing it won't change
            %            Output --
            %            Sig.test_signal Denoised signal
            %            Sig.Test_time After denoising the time vector may shift due to the convolution, the filter length of data get abandoned (normally less than 20)
            if nargin ==1
                CutoffFreq = [];
                FilterOrder = 20;
            end
            %if not specifiy the cutoff frequency, popout the stft plot
            if isempty(CutoffFreq)
                figure('Name',"Stft transfer of the signal")
                stft(theSignalAnalysis.Test_signal,theSignalAnalysis.Test_freq)
                fig = gcf;
                fig.Children(1).Label.String = 'Magnitude [dB]';
                ax = gca;
                X1 = ax.Children.XData;
                ax.Children.XData = X1*60;
                ax.XLim = [min(ax.Children.XData) Inf];
                ax.XLabel.String = 'Time [s]';
                ax.YLabel.String = 'Frequenct [Hz]';
                CutoffFreq = input('Set a Frequency High limits for the signal in Hz: ');
                if isempty(CutoffFreq)
                    CutoffFreq = 10;
                end
            end
            
%             [s_stft,f_stft,~] =...
%                 stft(theSignalAnalysis.Test_signal,theSignalAnalysis.Test_freq);
%             %istft
%             %set the frequency according to the STFT plot
%             stft_filter = find(abs(f_stft) > CutoffFreq);
%             % update 17-3 replace the other frequency by zero row
%             s_stft(stft_filter,:) = s_stft(stft_filter,:)*0;
%             [stft_signal,stft_t] = istft(s_stft,gather(theSignalAnalysis.Test_freq));
%             %Signal shift
%             stft_signal = real(stft_signal(15:end-15));
%             stft_t = stft_t(15:end-15);
%             
%             if ismember(Plotstates,theSignalAnalysis.BinaryY)
%                 figure('Name','The signal with stft denoised signal');
%                 t = tiledlayout(2,1);
%                 ax1 = nexttile;
%                 plot(theSignalAnalysis.Test_time,theSignalAnalysis.Test_signal);
%                 xlim([0 gather(max(theSignalAnalysis.Test_time))]);
%                 title('Original')
%                 ax2 = nexttile;
%                 plot(stft_t,stft_signal);
%                 xlim([0 gather(max(theSignalAnalysis.Test_time))]);
%                 title('Denoised STFT')
%                 linkaxes([ax1 ax2],'xy')
%                 ylabel(t, 'Current [pA]')
%                 xlabel(t, 'Time [s]')
%             end
%             
%             theSignalAnalysis.Test_signal = real(stft_signal);
%             theSignalAnalysis.Test_time = real(stft_t);
            % filter by the filter, hard to compute the scale factor
            Flag = 'Yes';
            N     = FilterOrder;  % Order
            while ismember(Flag,theSignalAnalysis.BinaryY)
                Fs = gather(theSignalAnalysis.Test_freq);  % Sampling Frequency
                
                Fpass = CutoffFreq ;  % Passband Frequency
                Fstop = CutoffFreq + 5;  % Stopband Frequency
                Wpass = 1;  % Passband Weight
                Wstop = 1;   % Stopband Weight
                
                % Calculate the coefficients using the FIRLS function.
                b  = firls(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop]);
                Hd = dfilt.dffir(b);
                
                fil_length = length(b);
                %try filtfilt with filter
                Sig_fil = filter(Hd,gather(theSignalAnalysis.Test_signal));
                %                 Sig_fil = filtfilt(Hd.Numerator,1,gather(theSignalAnalysis.Test_signal));
                Sig_fil = Sig_fil(fil_length:end);
                %resample the signal into the original length,effect by the
                %filter orders
                %                 Sig_fil = resample(Sig_fil,length(theSignalAnalysis.Test_signal),length(Sig_fil));
                scale_factor= mean([min(theSignalAnalysis.Test_signal)/min(Sig_fil),...
                    max(theSignalAnalysis.Test_signal)/max(Sig_fil)...
                    median(theSignalAnalysis.Test_signal)/median(Sig_fil)]);
                %                 RandIdx = randi(length(theSignalAnalysis.Test_signal),[1000,1]);
                %                 scale_factor = mean(theSignalAnalysis.Test_signal(RandIdx)./Sig_fil(RandIdx));
                figure('Name',sprintf('Cutoff frerquency:%d Hz,Order:%d \n',CutoffFreq,N));
                plot(theSignalAnalysis.Test_time(fil_length:end),Sig_fil*scale_factor,'DisplayName','Filtered Sig');
                hold on
                plot(theSignalAnalysis.Test_time,theSignalAnalysis.Test_signal,'DisplayName','Original')
                hold off;
                title(sprintf('Cutoff frerquency:%d Hz,Order:%d \n',CutoffFreq,N))
                xlabel("Time[s]")
                ylabel("Current [pA]")
                xlim([0 gather(max(theSignalAnalysis.Test_time))])
                legend show
                legend box off
                legend('Location', 'bestoutside')
                Flag = input('if change the order for the filter, input yes else no: ',"s");
                if ismember(Flag,theSignalAnalysis.BinaryY)
                    N = input('The order to change: ');
                end
            end
            theSignalAnalysis.Test_time = theSignalAnalysis.Test_time(fil_length:end);
            theSignalAnalysis.Test_signal = Sig_fil*scale_factor;
        end
        
        function [theSignalAnalysis] = BgSub(theSignalAnalysis,WinSize,Plotstates)
            %[theSignalAnalysis] = BgSub(theSignalAnalysis,WinSize,Plotstates)
            %offset the background by rloess
%             Input --
%             WinSize Windows size, by default set as 1000, lower the value, and the curve fits better the NIE spikes(undesired).
%             Plotstates whether or not showing the plot for the background subtraction, Output a subplot the first is the test_signal with a fitted trend line. The second one is the offset signal
%             Output--
%             Sig.test_signal_offset Offset signal data vector
%             Sig.test_signal_trendline Fitted trendline
            if nargin ==1
                WinSize = 1000;
                Plotstates = "N";
            end
            Trendline = smoothdata(gather(theSignalAnalysis.Test_signal),'rloess',WinSize);
            SigOff = theSignalAnalysis.Test_signal - Trendline;
            if ismember(Plotstates,theSignalAnalysis.BinaryY)
                figure;
                tiledlayout(2,2);
                ax1 = nexttile(1,[1,2]);
                plot(theSignalAnalysis.Test_time,theSignalAnalysis.Test_signal,'DisplayName','Original Signal')
                hold on;
                plot(theSignalAnalysis.Test_time,Trendline,'DisplayName','Trendline')
                hold off;
                ylabel( 'Current [pA]')
                xlabel( 'Time [s]')
                xlim([0 gather(max(theSignalAnalysis.Test_time))])
                title(sprintf("Trendline and Signal, Winsize:%d",WinSize))
                legend show
                legend('boxoff','Location','best')
                ax2 = nexttile(3,[1,2]);
                plot(theSignalAnalysis.Test_time,SigOff)
                ylabel( 'Current [pA]')
                xlabel( 'Time [s]')
                xlim([0 gather(max(theSignalAnalysis.Test_time))])
                title('Offset Signal')
                linkaxes([ax1 ax2],'x')
                ylabel( 'Current [pA]')
                xlabel( 'Time [s]')
            end
            theSignalAnalysis.Test_signal_trendline = Trendline;
            theSignalAnalysis.Test_signal_offset = SigOff;
        end
        
        function [theSignalAnalysis] = FlipFindPeak(theSignalAnalysis,BackgroundSigOffset,ReactionStates,Plotstates,PlotIndex)
            %[theSignalAnalysis] = FlipFindPeak(theSignalAnalysis,BackgroundSigOffset,ReactionStates,Plotstates,PlotIndex)
            %Flip find the peaks
            %
            %             Input --
            %             BackgroundSigOffset Blank data vector of the corresponding experiments
            %             ReactionStates define the oxidative or reductive NIE spikes, as 'Oxi' or 'Red'
            %             Plotstates Whether or not to show the plot e.g. 'Y'
            %             PlotIndex to plot found interval on the offset or original signal "Original" or "Offset"
            %             Output --
            %             Sig.Reactionstates same as input ReactionStates
            %             Sig.BackGroundHeightThr Height threshold calculated according to the Blank data
            %             Sig.FlipFindInte Found interval by height threshold in data index, first column left point, second column right points.
            %             Sig.PeaksFound three colum vector first colum as height, second the peak location in index, third the width unit as second.
            %
            

            if nargin ==1
            end
            theSignalAnalysis.Reactionstates = ReactionStates;
            %H_low is the 20 - 80% std of the blank Signal
            H_low = gather(2*std(BackgroundSigOffset(...
                ceil(length(BackgroundSigOffset)*0.2):ceil(length(BackgroundSigOffset)*0.8))));
            %             % https://doi.org/10.1016/j.coelec.2020.05.010
            %             lowestCurrent = gather(1e12 * 2100*1.60 * 1e-19 * theSignalAnalysis.Test_freq);
            theSignalAnalysis.BackGroundHeightThr = H_low;
            %find peaks
            switch ReactionStates
                case 'Red'
                    %shift up and down for offset signal
                    %                     theSignalAnalysis.Test_signal_offset = theSignalAnalysis.Test_signal_offset - std(theSignalAnalysis.Test_signal_offset);
                    [pks1, loc1, Pkswidth] = findpeaks(-gather(theSignalAnalysis.Test_signal_offset),'MinPeakHeight', H_low);
                    %                     [~, loc2, ~] = findpeaks(gather(theSignalAnalysis.Test_signal_offset));
                case 'Oxi'
                    %                     theSignalAnalysis.Test_signal_offset = theSignalAnalysis.Test_signal_offset + std(theSignalAnalysis.Test_signal_offset);
                    [pks1, loc1, Pkswidth] = findpeaks(gather(theSignalAnalysis.Test_signal_offset),'MinPeakHeight', H_low);
                    %                     [~, loc2, ~] = findpeaks(-gather(theSignalAnalysis.Test_signal_offset),'MinPeakHeight', 0);
            end
            SignSig = sign(theSignalAnalysis.Test_signal_offset);
            Idxnz = find(diff(SignSig)~=0);
            gap = abs(theSignalAnalysis.Test_signal_offset(Idxnz(1)) - theSignalAnalysis.Test_signal_offset(Idxnz(1) + 1));
            %             val = min(abs(theSignalAnalysis.Test_signal_offset));
            loc3 = find(abs(theSignalAnalysis.Test_signal_offset) < gap/2);
            %Debug
            %             figure;plot(theSignalAnalysis.Test_time,theSignalAnalysis.Test_signal_offset)
            %             hold on
            %             scatter(theSignalAnalysis.Test_time(loc1),theSignalAnalysis.Test_signal_offset(loc1),'b')
            %             hold on
            %             scatter(theSignalAnalysis.Test_time(loc2),theSignalAnalysis.Test_signal_offset(loc2),'r')
            %             hold on
            %             scatter(theSignalAnalysis.Test_time(loc3),theSignalAnalysis.Test_signal_offset(loc3),'g')
            loc2 = gather(loc3);
            %find the peaks in the flipped peaks
            LeftsideRawFlipInte = [1;loc2(1:end)];
            RightsideRawFlipInte = [loc2(1:end);length(theSignalAnalysis.Test_time)];
            InteFlipSelected = zeros(length(loc1),2);
            for i = 1:length(RightsideRawFlipInte)
                for j = 1:length(loc1)
                    if loc1(j) >= LeftsideRawFlipInte(i) && loc1(j) <= RightsideRawFlipInte(i)
                        InteFlipSelected(j,1) = LeftsideRawFlipInte(i);
                        InteFlipSelected(j,2) = RightsideRawFlipInte(i);
                    end
                end
            end
            %if serval peaks in one interval take the heighest, and remove
            %the repeated part
            PeaksProFound = [pks1, loc1, Pkswidth];
            [LPval,LPidx,~] =  unique(InteFlipSelected(:,1),'stable');
            [RPval,RPidx,~] =  unique(InteFlipSelected(:,2),'stable');
            %Repeated index on the right side
            RPrIdx = setdiff(1:size(PeaksProFound,1),RPidx);
            %pepeated index on the left side
            LPridx = setdiff(1:size(PeaksProFound,1),LPidx);
            
            %take the highest in the peak vector
            for i = 1:size(RPrIdx,2)
                Vecidx = find(InteFlipSelected(:,1) == InteFlipSelected(LPridx(i),1));
                Numr = sum(InteFlipSelected(:,1) == InteFlipSelected(LPridx(i),1));
                [~,Maxidx] = max([PeaksProFound(Vecidx,1)]);
                %except for the heighest peak rest become zero
                PeaksProFound(Vecidx(1) - 1 + setdiff(1:Numr,Maxidx),:) = zeros(Numr - 1,3);
            end
            PeaksProFound = [nonzeros(PeaksProFound(:,1)),nonzeros(PeaksProFound(:,2)),nonzeros(PeaksProFound(:,3))];
            %The interval take the longest
            if isequal(RPrIdx,LPridx)
                % if the repeated the interval is equal, take the vector of both side
                InteFlipSelected = [LPval,RPval];
            else
                % if the common on one side need to iterate each common
                % interval to find the longest
                warning('Two sides of interval are not even \n')
                for i = 1:size(LPrIdx,2)
                    Vecidx = InteFlipSelected(:,1) == InteFlipSelected(LPridx(i),1);
                    [Maxval,Maxidx] = max(InteFlipSelected(Vecidx ,2));
                    %need to be tested
                    InteFlipSelected(Vecidx(1) - 1 + Maxidx,:) = [InteFlipSelected(LPridx(i),1), Maxval];
                    InteFlipSelected(Vecidx(1) - setdiff(1:length(InteFlipSelected(Vecidx ,2)),Maxidx),:) = zeros(length(InteFlipSelected(Vecidx ,2)) - 1,3);
                end
                for i = 1:size(RPrIdx,2)
                    Vecidx = InteFlipSelected(:,2) == InteFlipSelected(RPridx(i),1);
                    [Maxval,Maxidx] = max(InteFlipSelected(Vecidx ,1));
                    InteFlipSelected(Vecidx(1) - 1 + Maxidx,:) = [Maxval,InteFlipSelected(RPridx(i),1)];
                    InteFlipSelected(Vecidx(1) - setdiff(1:length(InteFlipSelected(Vecidx ,2)),Maxidx),:) = zeros(length(InteFlipSelected(Vecidx ,2)) - 1,3);;
                end
                InteFlipSelected = [nonzeros(InteFlipSelected(:,1)),nonzeros(InteFlipSelected(:,2))];
            end
            theSignalAnalysis.FlipFindInte = InteFlipSelected;
            theSignalAnalysis.PeaksFound = PeaksProFound;
            %Plot section
            if ismember(Plotstates,theSignalAnalysis.BinaryY)
                PlotX = InteFlipSelected;
                s = sprintf('The signal with unmerged intervals');
                switch PlotIndex
                    case 'Original'
                        figure('Name','Marked all the interval in the original signal');
                        plot(theSignalAnalysis.Test_time,theSignalAnalysis.Test_signal)
                        hold on
                        scatter(theSignalAnalysis.Test_time(theSignalAnalysis.PeaksFound(:,2)),theSignalAnalysis.Test_signal(theSignalAnalysis.PeaksFound(:,2)))
                        hold on
                        for i = 1:length(PlotX)
                            hold on
                            plot(theSignalAnalysis.Test_time(PlotX(i,1):PlotX(i,2)),theSignalAnalysis.Test_signal(PlotX(i,1):PlotX(i,2)),'linewidth',3)
                        end
                        xlim([0 gather(max(theSignalAnalysis.Test_time))]);
                        title(sprintf(s +"," + "Oringial"));
                        ylabel('Current [pA]')
                        xlabel('Time [s]')
                    case 'Offset'
                        figure('Name','Marked all the interval in the offset signal');
                        plot(theSignalAnalysis.Test_time,theSignalAnalysis.Test_signal_offset)
                        hold on
                        scatter(theSignalAnalysis.Test_time(theSignalAnalysis.PeaksFound(:,2)),theSignalAnalysis.Test_signal_offset(theSignalAnalysis.PeaksFound(:,2)))
                        for i = 1:length(PlotX)
                            hold on
                            plot(theSignalAnalysis.Test_time(PlotX(i,1):PlotX(i,2)),theSignalAnalysis.Test_signal_offset(PlotX(i,1):PlotX(i,2)),'linewidth',3)
                            hold on
                        end
                        xlim([0 gather(max(theSignalAnalysis.Test_time))]);
                        title(sprintf(s +"," + "Offset"));
                        ylabel('Current [pA]')
                        xlabel('Time [s]')
                end
            end
        end
        
        function [theSignalAnalysis] = GeRawTrainSet(theSignalAnalysis,NumPlot,SilPlot)
            %[theSignalAnalysis] = GeRawTrainSet(theSignalAnalysis,'NumPlot')
            %Generate the raw template features training set and signal
            %if NumPlot is yes it will show the elbow method
            %if SilPlot is yes it will show the silhouette plot
            %Input--
            %NumPlot Elbow plot to decide the cluster numbers
            %SilPlot Silhouette plot to decide centroid number
            %Output--
            %Sig.RawTemClusterFeatures A table of clustering features, order same as Sig.FlipFindInte
            %Sig.RawTemClusterFeaturesNorm Z-score normalized cluster features
            %Sig.RawTemFeaturesMu Corresponding mean value for each feature
            %Sig.RawTemFeaturesSigma Corresponding standard deviation value for each feature
            InteToCluster = theSignalAnalysis.FlipFindInte;
            InteToCluster = [InteToCluster, ones(size(InteToCluster,1),1)];
            PeakIntergration = zeros(size(InteToCluster,1),1);
            for i = 1:size(InteToCluster,1)
                PeakIntergration(i) = polyarea(theSignalAnalysis.Test_time(InteToCluster(i,1):InteToCluster(i,2)),theSignalAnalysis.Test_signal_offset(InteToCluster(i,1):InteToCluster(i,2)));
            end
            %Calculation of the range mean median std of the signal in each
            %interval
            PeakbasicProperties = zeros(size(InteToCluster,1),4);
            for i = 1:size(InteToCluster,1)
                t_inte = theSignalAnalysis.Test_time(InteToCluster(i,1):InteToCluster(i,2));
                sig_inte = theSignalAnalysis.Test_signal_offset(InteToCluster(i,1):InteToCluster(i,2));
                %range
                PeakbasicProperties(i,1) = abs(min(t_inte)-max(t_inte));
                %mean inte
                PeakbasicProperties(i,2) = mean(sig_inte);
                % median range
                PeakbasicProperties(i,3) = median(sig_inte);
                % std inte
                PeakbasicProperties(i,4) = std(sig_inte);
            end
            %Calculated the interval properties
            IntePeaksProperties = zeros(size(InteToCluster,1),max(InteToCluster(:,3))*3 + 1);
            IntePeaksProperties(:,1) = InteToCluster(:,3);
            PeaksOutTimeWidth = theSignalAnalysis.PeaksFound;
            PeaksOutTimeWidth(:,3) =  PeaksOutTimeWidth(:,3)*(1/theSignalAnalysis.Test_freq);
            PeaksOutTimeWidth(:,4) = linspace(1,size(PeaksOutTimeWidth,1),size(PeaksOutTimeWidth,1));
            for i = 1:size(InteToCluster,1)
                PeaksIntheinterval = find(PeaksOutTimeWidth(:,4) == i);
                PeakInfo = zeros(1,size(PeaksIntheinterval,1)*3);
                %Peak Loc
                PeakInfo(1:3:end) = abs(PeaksOutTimeWidth(PeaksIntheinterval,2)-InteToCluster(i,1))/abs(InteToCluster(i,1)-InteToCluster(i,2));
                %Peak Height
                PeakInfo(2:3:end) =  PeaksOutTimeWidth(PeaksIntheinterval,1);
                %Peak Width
                PeakInfo(3:3:end) = PeaksOutTimeWidth(PeaksIntheinterval,3);
                %small peaks amounts in the interval? thinking to add
                IntePeaksProperties(i,2:size(PeaksIntheinterval,1)*3 + 1) = PeakInfo;
            end
            %change the InteFeatures into the table for future calling
            InteFeatures = [PeakIntergration,PeakbasicProperties,IntePeaksProperties(:,2:end)];
            InteFeaturesNorm = zeros(size(InteFeatures));
            muInteFeatures = zeros(size(InteFeatures,2),1);
            sigmaInteFeatures = zeros(size(InteFeatures,2),1);
            %Normlized the inte features
            for i = 1:size(InteFeatures,2)
                [InteFeaturesNorm(:,i),muInteFeatures(i),sigmaInteFeatures(i)] = SignalAnalysis.featureNormalize(InteFeatures(:,i));
            end
            %Signal data in each interval
            theSignalAnalysis.TrainingSignalOri =  zeros(size(InteToCluster,1),max(max((InteToCluster(:,2)-InteToCluster(:,1)))));
            for i = 1:size(InteToCluster,1)
                theSignalAnalysis.TrainingSignalOri = (SignalAnalysis.Autofill(theSignalAnalysis.Test_signal_offset(InteToCluster(i,1):InteToCluster(i,2)),theSignalAnalysis.TrainingSignalOri',i))';
            end
            theSignalAnalysis.RawTemFeaturesMu = muInteFeatures;
            theSignalAnalysis.RawTemFeaturesSigma = sigmaInteFeatures;
            theSignalAnalysis.RawTemClusterFeatures =  array2table(InteFeatures,...
                'VariableNames',{'Area','TimeDuration','Mean','Midean','Std','RelativePeakLoc','Peakheight','PeakWidth'});
            theSignalAnalysis.RawTemClusterFeaturesNorm = InteFeaturesNorm;
            %the plot for selecting the numbers for the clustering
            %15 is the test num for the elbow method and the silhouette
            %plot
            SignalAnalysis.KmeansClusterNumSelection(theSignalAnalysis,InteFeaturesNorm,15,NumPlot,SilPlot);
        end
        
        function [theSignalAnalysis] = KmeansGeRawSigTem(theSignalAnalysis,ClusterNum,Templot,Staplot)
            %[theSignalAnalysis] = KmeansGeRawSigTem(theSignalAnalysis,ClusterNum,Templot,Staplot)
            %Generation for the raw templates
            %Kmeans the Raw templates
            %Input--
            %ClusterNum Centroid numbers selected by the method in the above section, if leave it as blank select automatically by the matlab evalclusters functions, but it's suggested to choose according to the plot
            %Templot whether or not to show the plot with raw templates
            %Staplot whether or not to show the bar plot with each raw templates
            %Output--
            %Sig.NumTemplateRawUnRegu Raw templates numbers
            %Sig.TemplateRawUnRegu a cell with each cell of corresponding raw templates data
            %Sig.TemplateRawUnReguCounts two columns vector, first column count numbers for corresponding raw templates, second column for the count percentage of total spikes
            if nargin == 1 || isempty(ClusterNum)
                if length(theSignalAnalysis.RawTemClusterFeaturesNorm)>=15
                    ClusterNum = evalclusters(theSignalAnalysis.RawTemClusterFeaturesNorm,'kmeans','silhouette','KList',1:15);
                    ClusterNum = ClusterNum.OptimalK;
                else
                    num = size(theSignalAnalysis.RawTemClusterFeaturesNorm,1)-1;
                    ClusterNum = evalclusters(theSignalAnalysis.RawTemClusterFeaturesNorm,'kmeans','silhouette','KList',1:num);
                    ClusterNum = ClusterNum.OptimalK;
                end
            end
            if isempty(theSignalAnalysis.gpuStates)
                [SignalInteIdx,FeatureCentroidNorm,~,~] = kmeans(theSignalAnalysis.RawTemClusterFeaturesNorm,ClusterNum,'MaxIter',2000,'Replicates',2000);
            else
                options = statset('UseParallel',true);
                [SignalInteIdx,FeatureCentroidNorm,~,~] = kmeans(theSignalAnalysis.RawTemClusterFeaturesNorm,ClusterNum,'MaxIter',2000,'Replicates',2000,'Options',options);
            end
            %Reconstruct the norm centroids by the mu and sigma
            IntecenRenorm = zeros(size(FeatureCentroidNorm));
            for i = 1:size(FeatureCentroidNorm,2)
                IntecenRenorm(:,i) = SignalAnalysis.NormReconstruct(FeatureCentroidNorm(:,i),theSignalAnalysis.RawTemFeaturesMu(i),theSignalAnalysis.RawTemFeaturesSigma(i));
            end
            %Generate the templates using the Func RoughTemplateGe
            TemplateRawGe = SignalAnalysis.RoughTemplateGe(...
                theSignalAnalysis,theSignalAnalysis.TrainingSignalOri,theSignalAnalysis.RawTemClusterFeatures,ClusterNum,SignalInteIdx);
            
            StacountInte = zeros(ClusterNum,2);
            for i = 1:ClusterNum
                StacountInte(i,:) = nnz(SignalInteIdx == i);
            end
            StacountInte(:,2) = StacountInte(:,1)/sum(StacountInte(:,1));
            theSignalAnalysis.Intecolors = turbo(ClusterNum);
            theSignalAnalysis.TemplateRawUnRegu = TemplateRawGe;
            theSignalAnalysis.TemplateRawUnReguCounts = StacountInte;
            theSignalAnalysis.NumTemplateRawUnRegu = ClusterNum;
            %Plot section
            if ismember(Templot,theSignalAnalysis.BinaryY)
                figure;
                f = gobjects(ClusterNum,1);
                xmax = zeros(ClusterNum,1);
                ymax = zeros(ClusterNum,1);
                ymin = zeros(ClusterNum,1);
                for i = 1:ClusterNum
                    f(i) = subplot(ceil(ClusterNum/2),2,i);
                    Xinte = [0:(1/theSignalAnalysis.Test_freq):(size(TemplateRawGe{i},1)-1)*(1/theSignalAnalysis.Test_freq)]';
                    plot(Xinte,TemplateRawGe{i,:},'Color', theSignalAnalysis.Intecolors(i,:), 'Linewidth',2)
                    xmax(i) = max(Xinte);
                    ymax(i) = max(TemplateRawGe{i,:});
                    ymin(i) = min(TemplateRawGe{i,:});
                    title(sprintf('Clustered Interval %d, count Numbers %d',i,StacountInte(i,1)))
                end
                xlabel(f(1:ClusterNum),'Time [s]');
                ylabel(f(1:ClusterNum),'Current [pA]');
                ylim(f(1:ClusterNum),[min(ymin) max(ymax)]);
                xlim(f(1:ClusterNum),[0 max(xmax)]);
            end
            if ismember(Staplot,theSignalAnalysis.BinaryY)
                figure;
                label = 'Interval: ' + string(1:ClusterNum);
                p = bar(StacountInte(:,1),'FaceColor','flat');
                for i = 1:ClusterNum
                    p.CData(i,:) = theSignalAnalysis.Intecolors(i,:);
                    text(i,StacountInte(i,1),[num2str(StacountInte(i,2)*100,'%.3f') + "%"],'vert','bottom','horiz','center','Fontweight','bold','Fontsize',12);
                    box off
                end
                set(gca,'xticklabel',label)
                ylabel('Numbers of intervals')
            end
        end
        
        function [theSignalAnalysis] = RawtemplatesReguFunc(theSignalAnalysis,ManSelectedNum,Plotstates)
            %[theSignalAnalysis] = RawtemplatesReguFunc(theSignalAnalysis,ManSelectedNum,Plotstates)
            %Regulated the raw template, take from left side min to peak to
            %right side min
            %Removing the noisy templates and making raw templates fine
            %Input--
            %ManSelectedNum a row vector of selected numbers of raw templates to perform the fine process.if it leaves as blank it will only take the raw templates more than 15%
            %Plotstates Whether or not to show the fine templates plot e.g. 'y'
            %Output--
            %Sig.TemplatateRawRegu A cell with fine templates
            %Sig.TemplatateRawReguNum A vector with selected number
            if nargin == 1 || isempty(ManSelectedNum)
                ManSelectedNum = find(theSignalAnalysis.TemplateRawUnReguCounts(:,2)>0.15);
                ManSelectedNum = ManSelectedNum';
            end
            [Regulatedtemplates,RegulatedtemplatesNum] = SignalAnalysis.templatesRegu(...
                theSignalAnalysis,theSignalAnalysis.TemplateRawUnRegu,...
                theSignalAnalysis.TemplateRawUnReguCounts,ManSelectedNum);
            theSignalAnalysis.TemplatateRawRegu = Regulatedtemplates;
            theSignalAnalysis.TemplatateRawReguNum = RegulatedtemplatesNum;
            %Plot section
            if ismember(Plotstates,theSignalAnalysis.BinaryY)
                figure('Name','Regulated Raw templates');
                f = gobjects(length(RegulatedtemplatesNum),1);
                xmax = zeros(length(RegulatedtemplatesNum),1);
                ymax = zeros(length(RegulatedtemplatesNum),1);
                ymin = zeros(length(RegulatedtemplatesNum),1);
                for i = 1:length(RegulatedtemplatesNum)
                    f(i) = subplot(ceil(length(RegulatedtemplatesNum)/2),2,i);
                    Xinte = [0:(1/gather(theSignalAnalysis.Test_freq)):(size(Regulatedtemplates{i},1)-1)*(1/gather(theSignalAnalysis.Test_freq))]';
                    plot(Xinte,Regulatedtemplates{i,:}, 'Color', theSignalAnalysis.Intecolors(RegulatedtemplatesNum(i),:),'Linewidth',2)
                    xmax(i) = Xinte(end);
                    ymax(i) = max(Regulatedtemplates{i,:});
                    ymin(i) = min(Regulatedtemplates{i,:});
                    title(sprintf('Clustered Raw Template %d, count Numbers %d',RegulatedtemplatesNum(i),theSignalAnalysis.TemplateRawUnReguCounts(RegulatedtemplatesNum(i),1)))
                end
                xlabel(f(1:length(RegulatedtemplatesNum)),'Time [s]');
                ylabel(f(1:length(RegulatedtemplatesNum)),'Current [pA]');
                ylim(f(1:length(RegulatedtemplatesNum)),[gather(min(ymin)) gather(max(ymax))]);
                xlim(f(1:length(RegulatedtemplatesNum)),[0 gather(max(xmax))]);
            end
        end
        
        function [theSignalAnalysis] = Templatematching(theSignalAnalysis,Plotstates)
            %[theSignalAnalysis] = Templatematching(theSignalAnalysis,Plotstates)
            %Perform the template matching with the fine templates
            %[theSignalAnalysis] = Templatematching(theSignalAnalysis,Plotstates)
            %--Input
            %Plotstates Whether or not to show the similarity stem plots
            %-Output
            %Sig.SimilarityLag it's a cell. The first dimension is the
            %template's number, the first column is the cosine similarity
            %of each template, and the second is the corresponding interval for each similarity point in time. The left point defines the interval plus the template length as the right point.
            Similag = cell(size(theSignalAnalysis.TemplatateRawRegu,1),2);
            %the gpuarray is not proficent for the cross-correlation
            theSignalAnalysis.Test_signal_offset = gather(theSignalAnalysis.Test_signal_offset);
            theSignalAnalysis.Test_time = gather(theSignalAnalysis.Test_time);
            for i = 1:size(theSignalAnalysis.TemplatateRawRegu,1)
                Similag{i,1} = zeros(size(theSignalAnalysis.Test_signal_offset,1)-size(theSignalAnalysis.TemplatateRawRegu{i},1) + 1,1);
                Similag{i,2} = zeros(size(theSignalAnalysis.Test_signal_offset,1)-size(theSignalAnalysis.TemplatateRawRegu{i},1) + 1,2);
                % Normolized corr coef
                TemNorm = SignalAnalysis.featureNormalize(theSignalAnalysis.TemplatateRawRegu{i});
                % Normalized corr
                %                 TemNorm = theSignalAnalysis.TemplateSelected{i};
                for j = 1:(size(theSignalAnalysis.Test_signal_offset,1)-size(theSignalAnalysis.TemplatateRawRegu{i},1) + 1)
                    % Normolized corr coef
                    SigNorm = gather(SignalAnalysis.featureNormalize(theSignalAnalysis.Test_signal_offset(j:j+size(theSignalAnalysis.TemplatateRawRegu{i},1)-1)));
                    % Normalized corr
                    %                     SigNorm = theSignalAnalysis.Test_signal_offset(j:j+size(theSignalAnalysis.TemplateSelected{i},1)-1);
                    Similag{i,1}(j) = (SigNorm' * TemNorm)/...
                        (norm(SigNorm)*norm(TemNorm));
                    Similag{i,2}(j,1) = theSignalAnalysis.Test_time(j);
                    Similag{i,2}(j,2) = theSignalAnalysis.Test_time(j + size(theSignalAnalysis.TemplatateRawRegu{i},1) - 1);
                end
            end
            theSignalAnalysis.SimilarityLag = Similag;
            %Plot Section
            if ismember(Plotstates,theSignalAnalysis.BinaryY)
                figure;
                ymax = zeros(size(theSignalAnalysis.TemplatateRawRegu,1),1);
                ymin = zeros(size(theSignalAnalysis.TemplatateRawRegu,1),1);
                fsim = gobjects(size(theSignalAnalysis.TemplatateRawRegu,1),1);
                for i = 1:size(theSignalAnalysis.TemplatateRawRegu,1)
                    fsim(i) =  subplot(size(theSignalAnalysis.TemplatateRawRegu,1),1,i);
                    stem(Similag{i,2}(:,1),Similag{i,1},'MarkerSize',1,'Color',theSignalAnalysis.Intecolors(theSignalAnalysis.TemplatateRawReguNum(i),:))
                    title(sprintf('Template %d NormCrossCorrelationCoeff,Templates Counts %d',i,theSignalAnalysis.TemplateRawUnReguCounts(theSignalAnalysis.TemplatateRawReguNum(i),1)))
                    ymax(i) = max(Similag{i,1});
                    ymin(i) = min(Similag{i,1});
                end
                xlabel(fsim(1:size(theSignalAnalysis.TemplatateRawRegu,1)),'Time starting point [s]');
                ylabel(fsim(1:size(theSignalAnalysis.TemplatateRawRegu,1)),'Norm Similarity Coeff');
                ylim(fsim(1:size(theSignalAnalysis.TemplatateRawRegu,1)),[min(ymin) max(ymax)]);
                xlim(fsim(1:size(theSignalAnalysis.TemplatateRawRegu,1)),[0 Inf])
                linkaxes(fsim(1:size(theSignalAnalysis.TemplatateRawRegu,1)),'xy')
            end
            
        end
        function [theSignalAnalysis] = TemplatematchingFiltering(theSignalAnalysis,SimilarityLevel,StdFiltercoeff,HeightWidthratioCoeff,EachTemplateMatchCurve,TotalTemplateMatchCurve)
            %[theSignalAnalysis] = TemplatematchingFiltering(theSignalAnalysis,SimilarityLevel,StdFiltercoeff,HeightWidthratioCoeff,EachTemplateMatchCurve,TotalTemplateMatchCurev)
            %filter some matched signal according to the features of the
            %templates, and merge the multiply matches
            % create a cell, first col is the matched two sides by template
            % mtaching, second col is the filtered by std and ratio
            %             --Input
            %             SimilarityLevel Similarity threshold by default is set as 0.9
            %             StdFiltercoeff, computed the correlated template's standard deviation(Temstd) and compared the matched span's standard deviation with StdFiltercoeff * Temstd. If it is smaller, then filter out the span. Mainly filter extreme short spikes. by default set as 0.35
            %             HeightWidthratioCoeff, computed the ratio of the template's height and width ratio (TemRatio) and compared the matched spa's height-width ratio n with HeightWidthratioCoeff * TemRatio. If it is smaller, then filter out the span. Mainly filter extremely uneven spikes. by default as 0.35
            %             EachTemplateMatchCurve Whether or not to show the span matched by each template, color same as the generated templates.
            %             TotalTemplateMatchCurve Whether or not to show all merged spans of templates on the curve.
            %             --Output
            %             Sig.TemplatateMatchedInte A three-column matrix matched interval, first column left point, second column right point, third column peak point, they are all in indices of the sig.test_signal.
            %             Sig.AMInteSig The data in each interval, the shorter span, get fed with zero to reach the same length.
            %
            PeaksTim = cell(size(theSignalAnalysis.TemplatateRawRegu,1),2);
            for i = 1:size(theSignalAnalysis.TemplatateRawRegu,1)
                [~,PeaksIdx] = findpeaks(theSignalAnalysis.SimilarityLag{i,1},"MinPeakProminence",SimilarityLevel);
                PeaksTim{i,1} = theSignalAnalysis.SimilarityLag{i,2}(PeaksIdx,:);
                % second col of the cell is the template mayched interval
                % two sides index, and the max index
                PeaksTim{i,2} = zeros(length(PeaksTim{i,1}),3);
                Stdfilter = std((theSignalAnalysis.TemplatateRawRegu{i}));
                %Ratio filter, height divide the range
                switch theSignalAnalysis.Reactionstates
                    case 'Oxi'
                        [MaxTem,~] = max(theSignalAnalysis.TemplatateRawRegu{i});
                        %                         MinGapFilter = abs((theSignalAnalysis.TemplatateRawRegu{i}(end) - theSignalAnalysis.TemplatateRawRegu{i}(1))/MaxInte) ;
                        HeightWidthRatioFilter = abs(MaxTem - min(theSignalAnalysis.TemplatateRawRegu{i}))/(length(theSignalAnalysis.TemplatateRawRegu{i}) * (1/theSignalAnalysis.Test_freq));
                        %                         LRWidthRatioFilter = (MaxIdx - 1)/(length(theSignalAnalysis.TemplatateRawRegu{i}) - MaxIdx - 1);
                    case 'Red'
                        [MaxTem,~] = min(theSignalAnalysis.TemplatateRawRegu{i});
                        %                         MinGapFilter = abs((theSignalAnalysis.TemplatateRawRegu{i}(end) - theSignalAnalysis.TemplatateRawRegu{i}(1))/MaxInte) ;
                        HeightWidthRatioFilter = abs(max(theSignalAnalysis.TemplatateRawRegu{i}) - MaxTem )/(length(theSignalAnalysis.TemplatateRawRegu{i}) * (1/theSignalAnalysis.Test_freq));
                        %                         LRWidthRatioFilter = (MaxIdx - 1)/(length(theSignalAnalysis.TemplatateRawRegu{i}) - MaxIdx - 1);
                end
                %exam each interval
                for j = 1:size(PeaksTim{i,1},1)
                    switch theSignalAnalysis.Reactionstates
                        case 'Oxi'
                            %find the time number loc of the value
                            [~,PeaksTim{i,2}(j,1)] = min(abs(PeaksTim{i,1}(j,1) - theSignalAnalysis.Test_time));
                            [~,PeaksTim{i,2}(j,2)] = min(abs(PeaksTim{i,1}(j,2) - theSignalAnalysis.Test_time));
                            %find the max in the interval
                            [MaxInteVal,MaxInteIdx] = max(theSignalAnalysis.Test_signal_offset(PeaksTim{i,2}(j,1):PeaksTim{i,2}(j,2)));
                            %in case the max value is the last value of the
                            %interval
                            if MaxInteIdx + PeaksTim{i,2}(j,1)- 1 >= PeaksTim{i,2}(j,2)||MaxInteIdx + PeaksTim{i,2}(j,1)- 1<= PeaksTim{i,2}(j,1)
                                [pks,loc] = findpeaks((theSignalAnalysis.Test_signal_offset(PeaksTim{i,2}(j,1):PeaksTim{i,2}(j,2))));
                                [MaxInteVal,idxFindMax] = max(pks);
                                MaxInteIdx = loc(idxFindMax);
                            end
                            if isempty(MaxInteVal)
                                PeaksTim{i,2}(j,1) = 0;
                                PeaksTim{i,2}(j,2) = 0;
                                PeaksTim{i,2}(j,3) = 0;
                                continue
                                
                            end
                            [MinLeftVal,RealtiveLeftLoc] = min(theSignalAnalysis.Test_signal_offset(PeaksTim{i,2}(j,1):(PeaksTim{i,2}(j,1) + MaxInteIdx - 1)));
                            [MinRightVal,RealtiveRightLoc] = min(theSignalAnalysis.Test_signal_offset((PeaksTim{i,2}(j,1) + MaxInteIdx - 1):PeaksTim{i,2}(j,2)));
                            StdInte = std(theSignalAnalysis.Test_signal_offset((PeaksTim{i,2}(j,1)  + RealtiveLeftLoc -1):(PeaksTim{i,2}(j,1) + MaxInteIdx + RealtiveRightLoc - 1)));
                            HeightWidthRatioInte = (MaxInteVal - max([MinLeftVal,MinRightVal]))/((RealtiveRightLoc + MaxInteIdx  - RealtiveLeftLoc) * (1/theSignalAnalysis.Test_freq));
                        case 'Red'
                            %find the time number loc of the value
                            [~,PeaksTim{i,2}(j,1)] = min(abs(PeaksTim{i,1}(j,1) - theSignalAnalysis.Test_time));
                            [~,PeaksTim{i,2}(j,2)] = min(abs(PeaksTim{i,1}(j,2) - theSignalAnalysis.Test_time));
                            %find the max,in RED is min in the interval
                            [MaxInteVal,MaxInteIdx] = min(theSignalAnalysis.Test_signal_offset(PeaksTim{i,2}(j,1):PeaksTim{i,2}(j,2)));
                            %avoid the max is at the last in the interval
                            if MaxInteIdx + PeaksTim{i,2}(j,1)- 1 >= PeaksTim{i,2}(j,2)
                                [pks,loc] = findpeaks(-(theSignalAnalysis.Test_signal_offset(PeaksTim{i,2}(j,1):PeaksTim{i,2}(j,2))));
                                [MaxInteVal,idxFindMax] = max(abs(pks));
                                MaxInteIdx = loc(idxFindMax);
                            end
                            if isempty(MaxInteVal)
                                PeaksTim{i,2}(j,1) = 0;
                                PeaksTim{i,2}(j,2) = 0;
                                PeaksTim{i,2}(j,3) = 0;
                                continue
                            end
                            [MinLeftVal,RealtiveLeftLoc] = max(theSignalAnalysis.Test_signal_offset(PeaksTim{i,2}(j,1):(PeaksTim{i,2}(j,1) + MaxInteIdx - 1)));
                            [MinRightVal,RealtiveRightLoc] = max(theSignalAnalysis.Test_signal_offset((PeaksTim{i,2}(j,1) + MaxInteIdx - 1):PeaksTim{i,2}(j,2)));
                            StdInte = std(theSignalAnalysis.Test_signal_offset((PeaksTim{i,2}(j,1)  + RealtiveLeftLoc -1):(PeaksTim{i,2}(j,1) + MaxInteIdx + RealtiveRightLoc - 1)));
                            %                             MinGApInte = abs(MinLeftVal - MinRightVal)/abs(MaxInte);
                            HeightWidthRatioInte = abs(MaxInteVal - min([MinLeftVal,MinRightVal]))/((RealtiveRightLoc + MaxInteIdx  - RealtiveLeftLoc) * (1/theSignalAnalysis.Test_freq));
                            
                    end
                    %if the thresold exsist 
                    try
                        ThrCond =  (StdInte <= (Stdfilter *StdFiltercoeff))||(HeightWidthRatioInte <= (HeightWidthRatioFilter * HeightWidthratioCoeff));
                    catch
                        %                         warning('Judging conditions not fitted, Templates(i):%d,Interval(j):%d',i,j)
                        PeaksTim{i,2}(j,1) = 0;
                        PeaksTim{i,2}(j,2) = 0;
                        PeaksTim{i,2}(j,3) = 0;
                        continue
                    end
                    % if smaller the idx become zero otherwise remian
                    if ThrCond
                        PeaksTim{i,2}(j,1) = 0;
                        PeaksTim{i,2}(j,2) = 0;
                        PeaksTim{i,2}(j,3) = 0;
                    else
                        PeaksTim{i,2}(j,3) = PeaksTim{i,2}(j,1) + MaxInteIdx - 1;
                        % maybe -2
                        PeaksTim{i,2}(j,2) = PeaksTim{i,2}(j,1) + MaxInteIdx + RealtiveRightLoc - 1;
                        PeaksTim{i,2}(j,1) = PeaksTim{i,2}(j,1) + RealtiveLeftLoc -1;
                    end
                    %last run remove the zero
                    if j == length(PeaksTim{i,1})
                        PeaksTim{i,2} = [nonzeros(PeaksTim{i,2}(:,1)),nonzeros(PeaksTim{i,2}(:,2)),nonzeros(PeaksTim{i,2}(:,3))];
                    end
                end
            end
            [~,InteMatchOrder] = sort(cellfun(@numel,PeaksTim(:,2)),'descend');
            
            for i = 1:size(theSignalAnalysis.TemplatateRawRegu,1)
                %created a exam vector, i =1; is made by other intervals selected by other templates
                %else
                %made by the vector in the lower order removed the row has same max with
                % finalVec
                % examvec to be scaned by testVec, the InteVec is mergedvec
                % the longest vec
                examVector = [];
                for j = i:size(theSignalAnalysis.TemplatateRawRegu,1)
                    if i == 1
                        examVector = [examVector;PeaksTim{InteMatchOrder(j),2}];
                    else
                        %intersectionidx of the j interval
                        [~,DiffFinMaxInteIdx]  = setdiff(PeaksTim{InteMatchOrder(j),2}(:,3),FinalVec(:,3));
                        examVector = [examVector;PeaksTim{InteMatchOrder(j),2}(DiffFinMaxInteIdx,:)];
                    end
                    %                     PeaksTim{InteMatchOrder(i),2}(:,3)
                end
                if i == 1
                    testVec = PeaksTim{InteMatchOrder(i),2};
                    InteVec = zeros(size(PeaksTim{InteMatchOrder(i),2}));
                    FinalVec = [];
                else
                    %find different elements in the array
                    [~,DiffTestFinalIdx]  = setdiff(PeaksTim{InteMatchOrder(i),2}(:,3),FinalVec(:,3));
                    testVec = PeaksTim{InteMatchOrder(i),2}(DiffTestFinalIdx,:);
                    InteVec = zeros(size(testVec));
                end
                %iteration for every intervals in found interval i
                for k = 1:size(testVec,1)
                    flag = 0;
                    for m = 1:size(examVector,1)
                        if ~ismember(testVec(k,3),examVector(:,3))
                            InteVec(k,:) = examVector(m,:);
                            continue
                        else
                            if testVec(k,3) == examVector(m,3) && flag == 0
                                InteVec(k,1) = min([testVec(k,1),examVector(m,1)]);
                                InteVec(k,2) = max([testVec(k,2),examVector(m,2)]);
                                InteVec(k,3) = testVec(k,3);
                                flag = flag + 1;
                            elseif testVec(k,3) == examVector(m,3) && flag ~= 0
                                InteVec(k,1) =  min([InteVec(k,1),examVector(m,1)]);
                                InteVec(k,2) =  max([InteVec(k,2),examVector(m,2)]);
                                flag = flag + 1;
                            else
                                continue
                            end
                            
                        end
                    end
                end
                FinalVec = [FinalVec;InteVec];
            end
            %FinalVec needed to to reselected, elimiate the left than
            %higher and the small interval in the big one itself
            [~,FinInteIdx] = sort(FinalVec(:,3));
            FinalVec = FinalVec(FinInteIdx,:);
            %examing two cases when there are overlaped two intervals which
            %there don't share the same max point
            i = 2;
            while i < size(FinalVec,1)
                if FinalVec(i,1) > FinalVec(i - 1,1) && FinalVec(i,1) < FinalVec(i - 1,2)
                    Intersect = intersect([FinalVec(i - 1,1):FinalVec(i - 1,2)],[FinalVec(i,1):1:FinalVec(i,2)]);
                    if FinalVec(i - 1,2) < FinalVec(i,3)
                        Midpoint = max(Intersect);
                    else
                        Midpoint = min(Intersect);
                    end
                    FinalVec(i - 1,2) = Midpoint;
                    FinalVec(i,1) = Midpoint;
                    i = i + 1;
                else
                    i = i + 1;
                end
            end
            FinalVec = [nonzeros(FinalVec(:,1)), nonzeros(FinalVec(:,2)),nonzeros(FinalVec(:,3))];
            theSignalAnalysis.TemplatateMatchedInte = FinalVec;
            %Take the signal in the matched interval
            theSignalAnalysis.AMInteSig =  zeros(size(FinalVec,1),max(max((FinalVec(:,2)-FinalVec(:,1)))));
            for i = 1:size(FinalVec,1)
                theSignalAnalysis.AMInteSig = (SignalAnalysis.Autofill(theSignalAnalysis.Test_signal_offset(FinalVec(i,1):FinalVec(i,2)),theSignalAnalysis.AMInteSig',i))';
            end
            %Plot Section
            if ismember(EachTemplateMatchCurve,theSignalAnalysis.BinaryY)
                figure;
                plot(theSignalAnalysis.Test_time,theSignalAnalysis.Test_signal_offset);
                xlim([min(theSignalAnalysis.Test_time) max(theSignalAnalysis.Test_time)])
                title('Offset signal with the matched interval by each template')
                xlabel('Time [s]')
                ylabel('Current [pA]')
                hold on
                %add another loop to iterate the different templates
                for i = 1:size(theSignalAnalysis.TemplatateRawRegu,1)
                    for j = 1:size(PeaksTim{i,2},1)
                        
                        plot(theSignalAnalysis.Test_time...
                            (PeaksTim{i,2}(j,1):PeaksTim{i,2}(j,2)),theSignalAnalysis.Test_signal_offset(PeaksTim{i,2}(j,1):PeaksTim{i,2}(j,2)),...
                            'linewidth',3,'Color',theSignalAnalysis.Intecolors(theSignalAnalysis.TemplatateRawReguNum(i),:))
                        
                    end
                end
            end
            
            if ismember(TotalTemplateMatchCurve,theSignalAnalysis.BinaryY)
                figure;
                plot(theSignalAnalysis.Test_time,theSignalAnalysis.Test_signal_offset);
                xlim([min(theSignalAnalysis.Test_time) max(theSignalAnalysis.Test_time)])
                title('Original signal with the all merged template matched interval')
                xlabel('Time [s]')
                ylabel('Current [pA]')
                hold on
                for i = 1:size(theSignalAnalysis.TemplatateMatchedInte,1)
                    plot(theSignalAnalysis.Test_time(theSignalAnalysis.TemplatateMatchedInte(i,1):theSignalAnalysis.TemplatateMatchedInte(i,2))...
                        ,theSignalAnalysis.Test_signal_offset(theSignalAnalysis.TemplatateMatchedInte(i,1):theSignalAnalysis.TemplatateMatchedInte(i,2)),'LineWidth',2)
                    hold on
                end
            end
        end
        function [theSignalAnalysis] = GeAMTrainSet(theSignalAnalysis,HeightFilter,NumPlot,SilPlot)
            %[theSignalAnalysis] = GeAMTrainSet(theSignalAnalysis)
            %Genearte the traning set for after the template matching
            %SlopeL,SlopeR,Range,RaltiveLoc
            DataOffset = zeros(size(theSignalAnalysis.TemplatateMatchedInte,1),5);
            % Height Area caluclate on the oringinal curve
            DataOri  = zeros(size(theSignalAnalysis.TemplatateMatchedInte,1),3);
            for i = 1:size(theSignalAnalysis.TemplatateMatchedInte,1)
                %SlopeL
                DataOffset(i,1) = theSignalAnalysis.Test_signal_offset(theSignalAnalysis.TemplatateMatchedInte(i,3)) - theSignalAnalysis.Test_signal_offset(theSignalAnalysis.TemplatateMatchedInte(i,1))/...
                    ((theSignalAnalysis.TemplatateMatchedInte(i,3) - theSignalAnalysis.TemplatateMatchedInte(i,1)) * (1/theSignalAnalysis.Test_freq));
                %SloprR
                DataOffset(i,2) = theSignalAnalysis.Test_signal_offset(theSignalAnalysis.TemplatateMatchedInte(i,2)) - theSignalAnalysis.Test_signal_offset(theSignalAnalysis.TemplatateMatchedInte(i,3))/...
                    ((theSignalAnalysis.TemplatateMatchedInte(i,2) - theSignalAnalysis.TemplatateMatchedInte(i,3)) * (1/theSignalAnalysis.Test_freq));
                %Range
                DataOffset(i,3) = (theSignalAnalysis.TemplatateMatchedInte(i,2) - theSignalAnalysis.TemplatateMatchedInte(i,1)) * (1/theSignalAnalysis.Test_freq);
                %Peak RelativeLoc
                DataOffset(i,4) = (theSignalAnalysis.TemplatateMatchedInte(i,3) - theSignalAnalysis.TemplatateMatchedInte(i,1))/(theSignalAnalysis.TemplatateMatchedInte(i,2) - theSignalAnalysis.TemplatateMatchedInte(i,1));
                %Height
                DataOffset(i,5) = abs(theSignalAnalysis.Test_signal_offset(theSignalAnalysis.TemplatateMatchedInte(i,3)));
                %the Prominence and area are calultaed on the oriignial plot
                switch theSignalAnalysis.Reactionstates
                    case 'Oxi'
                        DataOri(i,1) = theSignalAnalysis.Test_signal(theSignalAnalysis.TemplatateMatchedInte(i,3)) -...
                            min([theSignalAnalysis.Test_signal(theSignalAnalysis.TemplatateMatchedInte(i,1)) theSignalAnalysis.Test_signal(theSignalAnalysis.TemplatateMatchedInte(i,2))]) ;
                    case 'Red'
                        DataOri(i,1) = max([theSignalAnalysis.Test_signal(theSignalAnalysis.TemplatateMatchedInte(i,1)) theSignalAnalysis.Test_signal(theSignalAnalysis.TemplatateMatchedInte(i,2))]) - ...
                            theSignalAnalysis.Test_signal(theSignalAnalysis.TemplatateMatchedInte(i,3));
                end
                %area
                DataOri(i,2) = polyarea(theSignalAnalysis.Test_time(theSignalAnalysis.TemplatateMatchedInte(i,1):theSignalAnalysis.TemplatateMatchedInte(i,2)),...
                    theSignalAnalysis.Test_signal(theSignalAnalysis.TemplatateMatchedInte(i,1):theSignalAnalysis.TemplatateMatchedInte(i,2)));
                %Prominence range ratio
                DataOri(i,3) = DataOri(i,1)/DataOffset(i,3);
            end
            %Editing
            AFTInteProperities = [DataOffset,DataOri];
            if ismember(HeightFilter,theSignalAnalysis.BinaryY)
                %for the signal which is noisy or low signal to noise ratio
                theSignalAnalysis.TemplatateMatchedInte = theSignalAnalysis.TemplatateMatchedInte(AFTInteProperities(:,5)> theSignalAnalysis.BackGroundHeightThr,:); 
                AFTInteProperities = AFTInteProperities(AFTInteProperities(:,5)> theSignalAnalysis.BackGroundHeightThr,:);
            end
            theSignalAnalysis.AMTraningSet = array2table(AFTInteProperities,...
                'VariableNames',{'SlopeL','SlopeR','TimeDuration','RelativePeakLoc','Height','Prominence','Area','ProminenceRangeRatio'});
            Pmu = zeros(1,size(AFTInteProperities,2));
            Psigma = zeros(1,size(AFTInteProperities,2));
            AFTInteProperitiesNorm = zeros(size(AFTInteProperities));
            for i = 1:size(AFTInteProperities,2)
                [AFTInteProperitiesNorm(:,i),Pmu(i),Psigma(i)] = SignalAnalysis.featureNormalize(AFTInteProperities(:,i));
            end
            theSignalAnalysis.AMTraningSetNorm = AFTInteProperitiesNorm;
            theSignalAnalysis.AMTemPSigma = Psigma;
            theSignalAnalysis.AMTemPMu = Pmu;
            if size(AFTInteProperities,1) ~= 1
                SignalAnalysis.KmeansClusterNumSelection(theSignalAnalysis,AFTInteProperitiesNorm,15,NumPlot,SilPlot);
            else
                fprintf("Only one peak signal\n")
            end
        end
        
        function [theSignalAnalysis] = KmeansGeAMSigTem(theSignalAnalysis,ClusterNum,MarkPlot,TotalPlot)
            %[theSignalAnalysis] = KmeansGeAMSigTem(theSignalAnalysis,ClusterNum,MarkPlot,TemPlot,TotalPlot)
            %Clustering the intervals by template matching
            %Generating the after matching templates
            if size(theSignalAnalysis.AMTraningSetNorm,1) ~= 1
                [TemPidx,TemPCen,~,~] = kmeans(theSignalAnalysis.AMTraningSetNorm,ClusterNum,'MaxIter',1000,'Replicates',500);
            else
                TemPidx = 1; TemPCen = theSignalAnalysis.AMTraningSetNorm;
            end
            IntecenRenorm = zeros(size(TemPCen));
            for i = 1:size(TemPCen,2)
                IntecenRenorm(:,i) = SignalAnalysis.NormReconstruct(...
                    TemPCen(:,i), theSignalAnalysis.AMTemPMu(i),theSignalAnalysis.AMTemPSigma(i));
            end
            theSignalAnalysis.AMTemCen = IntecenRenorm;
            theSignalAnalysis.AMTemidx = TemPidx;
            TemplateRawGe = SignalAnalysis.RoughTemplateGe(theSignalAnalysis...
                ,theSignalAnalysis.AMInteSig,theSignalAnalysis.AMTraningSet,...
                ClusterNum,TemPidx);
            SelectNum = 1:1:ClusterNum;
            TemStacountInte = zeros(ClusterNum,2);
            for i = 1:ClusterNum
                TemStacountInte(i,1) = nnz(TemPidx == i);
            end
            TemStacountInte(:,2) = TemStacountInte(:,1)/sum(TemStacountInte(:,1));
            [Regulatedtemplates,~] = ...
                SignalAnalysis.templatesRegu(...
                theSignalAnalysis,TemplateRawGe,TemStacountInte,SelectNum);
            theSignalAnalysis.AMTemplates = Regulatedtemplates;
            theSignalAnalysis.AMIntecolors = turbo(ClusterNum);
            %PlotSection
            if ismember(MarkPlot,theSignalAnalysis.BinaryY)
                Inte = theSignalAnalysis.TemplatateMatchedInte;
                yplot = theSignalAnalysis.Test_signal;
                sg = 'Original plot';
                figure;
                plot(theSignalAnalysis.Test_time,yplot,'DisplayName',sg)
                hold on
                % all the lines in one type are ploted in one row
                pltinte = gobjects(ClusterNum,max(TemStacountInte(1,:)));
                for j = 1:ClusterNum
                    NumPLotsInte = size(Inte(TemPidx == j,1),1);
                    InteSidesPlot = zeros(NumPLotsInte,2);
                    InteSidesPlot(:,1) = Inte(TemPidx == j,1);
                    InteSidesPlot(:,2) = Inte(TemPidx == j,2);
                    sginte = ['Interval type: ',num2str(j)];
                    for k = 1:NumPLotsInte
                        pltinte(j,k) =  plot(theSignalAnalysis.Test_time(InteSidesPlot(k,1):InteSidesPlot(k,2)),yplot(InteSidesPlot(k,1):InteSidesPlot(k,2)),'linewidth',3,'Color',theSignalAnalysis.AMIntecolors(j,:),"DisplayName",sginte);
                        hold on
                    end
                end
                title('Signal with the clustered peaks with clustered interval')
                xlabel('Time [s]')
                ylabel('Current [pA]')
                xlim([0 max(theSignalAnalysis.Test_time)]);
                legend([pltinte(:,1)],'box','off','Location','northeastoutside')
                hold off
            end
            if ismember(TotalPlot,theSignalAnalysis.BinaryY)
                figure;
                tiledlayout( 4 + 2*ClusterNum , 8,"Padding","loose")
                %fisrt one is the sum plot
                nexttile(1,[4 2]);
                label = 'Type: ' + string(1:ClusterNum);
                p = bar(TemStacountInte(:,1),'FaceColor','flat');
                for i = 1:ClusterNum
                    p.CData(i,:) = theSignalAnalysis.AMIntecolors(i,:);
                    text(i,TemStacountInte(i,1),[num2str(TemStacountInte(i,2)*100,'%.3f') + "%"],'vert','bottom','horiz','center','FontSize',9);
                    box off
                end
                set(gca,'xticklabel',label)
                ylabel('Counts')
                %Area distribution
                nexttile(3,[4 2]);
                histogram(theSignalAnalysis.AMTraningSet.Area,ceil(sqrt(size(theSignalAnalysis.AMTraningSet,1))))
                xlabel('Area[pC]')
                ylabel('Counts')
                %Height
                nexttile(5,[4 2]);
                histogram(theSignalAnalysis.AMTraningSet.Height,ceil(sqrt(size(theSignalAnalysis.AMTraningSet,1))))
                ylabel('Counts')
                xlabel('Height[pA]')
                %range
                nexttile(7,[4 2]);
                histogram(theSignalAnalysis.AMTraningSet.TimeDuration,ceil(sqrt(size(theSignalAnalysis.AMTraningSet,1))))
                ylabel('Counts')
                xlabel('TimeDuration[s]')
                % Avgtemplates
                xmax = zeros(ClusterNum,1);
                ymax = zeros(ClusterNum,1);
                ymin = zeros(ClusterNum,1);
                for i = 1:ClusterNum
                    axTem(i) = nexttile(i*16 + 17,[2 2]);
                    Xinte = [0:(1/theSignalAnalysis.Test_freq):(size(Regulatedtemplates{i},1)-1)*(1/theSignalAnalysis.Test_freq)]';
                    plot(Xinte,Regulatedtemplates{i,:}, 'Color', theSignalAnalysis.AMIntecolors(i,:),'Linewidth',2)
                    xmax(i) = max(Xinte);
                    ymax(i) = max(Regulatedtemplates{i,:});
                    ymin(i) = min(Regulatedtemplates{i,:});
                    title(sprintf('Templates %d, count Numbers %d',i,TemStacountInte(i,1)))
                    ylim([ymin(i) ymax(i)* 1.05])
                end
                xlim(axTem,[0 max(xmax)])
                %                 ylim(axTem,[min(ymin) 1.05 * max(ymax)])
                xlabel(axTem,'Time [s]')
                ylabel(axTem,'Current [pA]')
                %Area Each Tem
                for i = 1:ClusterNum
                    axArea(i) = nexttile(i*16 + 19,[2 2]);
                    histogram(theSignalAnalysis.AMTraningSet.Area(TemPidx == i),ceil(sqrt( nnz(TemPidx == i))),'FaceColor',theSignalAnalysis.AMIntecolors(i,:))
                    title(sprintf('Area Distribution Type %d',i))
                end
                xlim(axArea,[0 max(theSignalAnalysis.AMTraningSet.Area)])
                xlabel(axArea,'Area[pC]')
                ylabel(axArea,'Counts')
                %Height Each Tem
                for i = 1:ClusterNum
                    axHeight(i) = nexttile(i*16 + 21,[2 2]);
                    histogram(theSignalAnalysis.AMTraningSet.Height(TemPidx == i),ceil(sqrt( nnz(TemPidx == i))),'FaceColor',theSignalAnalysis.AMIntecolors(i,:))
                    title(sprintf('Height Distribution Type %d',i))
                end
                xlim(axHeight,[0 max(theSignalAnalysis.AMTraningSet.Height)])
                xlabel(axHeight,'Height[pA]')
                ylabel(axHeight,'Counts')
                %Range
                for i = 1:ClusterNum
                    axRange(i) = nexttile(i*16 + 23,[2 2]);
                    histogram(theSignalAnalysis.AMTraningSet.TimeDuration(TemPidx == i),ceil(sqrt( nnz(TemPidx == i))),'FaceColor',theSignalAnalysis.AMIntecolors(i,:))
                    title(sprintf('Duration Distribution Type %d',i))
                end
                xlim(axRange,[0 max(theSignalAnalysis.AMTraningSet.TimeDuration)])
                xlabel(axRange,'TimeDuration[s]')
                ylabel(axRange,'Counts')
                
            end
        end
%         function [theSignalAnalysis] = ConGaussian(theSignalAnalysis,PlotState)
%             %[theSignalAnalysis] = ConGaussian(theSignalAnalysis,PlotState)
%             %Compoute the data of the original set up by guassain
%             %distribuion, will plot the comoarasion of the data computed
%             %assuming the peak is gaussain 
%             DGaussain = zeros(size(theSignalAnalysis.RawTemClusterFeatures,1),6);
%             DGaussain(:,1) = theSignalAnalysis.RawTemClusterFeatures.Peakheight;
%             %std computed by the width 
%             DGaussain(:,2) = theSignalAnalysis.RawTemClusterFeatures.PeakWidth/2.355;
%             %area H*std*0.3989
%             DGaussain(:,3) = (DGaussain(:,1) .* DGaussain(:,2))/0.3989;
%             % peaks location
%             DGaussain(:,5) = theSignalAnalysis.PeaksFound(:,2);
%             %left side peak point -3*sigma
%             DGaussain(:,4) = DGaussain(:,5) - round(3 * DGaussain(:,2)*gather(theSignalAnalysis.Test_freq));
%             DGaussain(DGaussain(:,4) <= 0,4) = 1;
%             %right point peak point +3*sigma
%             DGaussain(:,6) = DGaussain(:,5) + round(3 * DGaussain(:,2)*gather(theSignalAnalysis.Test_freq));
%             %             pks1(i,:) * (gaussmf(time_simInt,[Pkswidth(i)*(1/thePeakFind.Test_freq)/2.355 thePeakFind.Test_Time(loc1(i))]));
%             Gausig = cell(size(DGaussain,1),2);
%             for i = 1:size(Gausig,1)
%                 Gausig{i,1} = [theSignalAnalysis.Test_time(DGaussain(i,4)):gather(1/theSignalAnalysis.Test_freq):theSignalAnalysis.Test_time(DGaussain(i,6))];
%                 if theSignalAnalysis.Reactionstates == "Oxi"
%                     Gausig{i,2} = DGaussain(i,1) * gaussmf(Gausig{i,1},[DGaussain(i,2),theSignalAnalysis.Test_time(DGaussain(i,5))]);
%                 else
%                     Gausig{i,2} = -DGaussain(i,1) * gaussmf(Gausig{i,1},[DGaussain(i,2),theSignalAnalysis.Test_time(DGaussain(i,5))]);
%                 end
%             end
%             theSignalAnalysis.DataGaussain = DGaussain;
%             theSignalAnalysis.SigGaussain = Gausig;
%             %find the shared value
% %             [~,idxDG] = intersect(DGaussain(:,5),theSignalAnalysis.TemplatateMatchedInte(:,3),'stable');
% %             [~,idxTM] = intersect(theSignalAnalysis.TemplatateMatchedInte(:,3),DGaussain(:,5),'stable');
% %             %Compare cell, with area and range for guassain is set as 6 time std
% %             DCompare = cell(1,2);
% %             %area
% %             DCompare{1} = zeros(size(idxDG,1),3);
% %             DCompare{1}(:,1) = DGaussain(idxDG,3);
% %             DCompare{1}(:,2) = theSignalAnalysis.AMTraningSet.Area(idxTM);
% %             DCompare{1}(:,3) = DCompare{1}(:,1) - DCompare{1}(:,2);
% %             %histogram(DCompare{1}(:,3))
% %             %mean(DCompare{1}(:,3));
% %             %std()
% %             %time range
% %             DCompare{2} = zeros(size(idxDG,1),3);
% %             DCompare{2}(:,1) = DGaussain(idxDG,2)*6;
% %             DCompare{2}(:,2) = theSignalAnalysis.AMTraningSet.TimeDuration(idxTM);
% %             % in ms
% %             DCompare{2}(:,3) = (DCompare{2}(:,1) - DCompare{2}(:,2))*1e3;
%             %mean(DCompare{2}(:,3));
%             %NUMBERS DIFFERENCE
%             %size(theSignalAnalysis.AMTraningSet,1)
%             %size(DGaussain,1)
%             %num the old way in fact more?
%             %mean(size(DGaussain,1) - size(theSignalAnalysis.AMTraningSet,1)
%             if ismember(PlotState,theSignalAnalysis.BinaryY)
%                 figure;
%                 plot(theSignalAnalysis.Test_time,theSignalAnalysis.Test_signal_offset);
%                 xlim([min(theSignalAnalysis.Test_time) max(theSignalAnalysis.Test_time)])
%                 title('Offset signal with peaks found by template matching and height thresold')
%                 xlabel('Time [s]')
%                 ylabel('Current [pA]')
%                 hold on
%                 ax1 = matlab.graphics.chart.primitive.Line.empty(size(theSignalAnalysis.TemplatateMatchedInte,1),0);
%                 for i = 1:size(theSignalAnalysis.TemplatateMatchedInte,1)
%                     ax1(i) = plot(theSignalAnalysis.Test_time(theSignalAnalysis.TemplatateMatchedInte(i,1):theSignalAnalysis.TemplatateMatchedInte(i,2))...
%                         ,theSignalAnalysis.Test_signal_offset(theSignalAnalysis.TemplatateMatchedInte(i,1):theSignalAnalysis.TemplatateMatchedInte(i,2)),...
%                         'LineWidth',2,'Color',"black",'LineWidth',1.5,'DisplayName','Signal span found by Template Matching');
%                     hold on
%                 end
%                 axsc1 = scatter(theSignalAnalysis.Test_time(theSignalAnalysis.TemplatateMatchedInte(:,3)),theSignalAnalysis.Test_signal_offset(theSignalAnalysis.TemplatateMatchedInte(:,3)),'+','black','LineWidth',1.5,'SizeData',65,'DisplayName','Peaks found by template matching');
%                 hold on
%                 ax2 = matlab.graphics.chart.primitive.Line.empty(size(theSignalAnalysis.TemplatateMatchedInte,1),0);
%                 for i = 1:size(Gausig,1)
%                     ax2(i) = plot(Gausig{i,1},Gausig{i,2},'Color',"red",'LineWidth',2,'DisplayName','Assumed Gaussain Distribution');
%                     hold on
%                 end
%                 axsc2 = scatter(theSignalAnalysis.Test_time(DGaussain(:,5)),theSignalAnalysis.Test_signal_offset(DGaussain(:,5)),'red','SizeData',72,'LineWidth',1.5,'DisplayName','Peaks found by Height threshold');
%                 legend([axsc1,axsc2,ax1(1),ax2(1)],'Box','off') 
%             end
%             
%         end
        function [theSignalAnalysis] = Rematching(theSignalAnalysis,TimesRematching,MarkPlot,TotalPlot)
            % Remacthching the signal by replacing TemplatateRawRegu by AMTemplates
            N = length(theSignalAnalysis.AMTemplates);
            for i = 1:TimesRematching
                theSignalAnalysis.TemplateRawUnRegu = theSignalAnalysis.AMTemplates;
                theSignalAnalysis = theSignalAnalysis.Templatematching('N');
                theSignalAnalysis = theSignalAnalysis.TemplatematchingFiltering(0.95,0.35,0.35,'N','N');
                theSignalAnalysis = theSignalAnalysis.GeAMTrainSet('N','N',"N");
                if i == TimesRematching
                    theSignalAnalysis = theSignalAnalysis.KmeansGeAMSigTem(N,MarkPlot,TotalPlot);
                else
                    theSignalAnalysis = theSignalAnalysis.KmeansGeAMSigTem(N,"N","N");
                end
            end

        end
    end
    
    
    
    
    methods(Static)
        
        function MatrixFilled = Autofill(vector,matrix,j)
            %   matrix = Autofill(vector,matrix,j)
            %   fill the vector into a matrix at j coloum if length different
            %   automatically fill 0 at the smaller length part
            
            MatrixFilled_T = vector;
            MatrixFilled = matrix;
            if size(MatrixFilled,1) == length(MatrixFilled_T)
                MatrixFilled(:,j) = MatrixFilled_T;
            elseif size(MatrixFilled,1) > length(MatrixFilled_T)
                MatrixFilled_T = [MatrixFilled_T;zeros(size(MatrixFilled,1)-length(MatrixFilled_T),1)];
                MatrixFilled(:,j) = MatrixFilled_T;
            elseif size(MatrixFilled,1) < length(MatrixFilled_T)
                if j >= size(MatrixFilled,2)
                    MatrixFilled =  [MatrixFilled zeros(size(MatrixFilled,1),j-size(MatrixFilled,2));zeros(length(MatrixFilled_T)-size(MatrixFilled,1),j)];
                    MatrixFilled(:,j) = MatrixFilled_T;
                else
                    MatrixFilled =  [MatrixFilled; zeros(length(MatrixFilled_T)-size(MatrixFilled,1),size(MatrixFilled,2))];
                    MatrixFilled(:,j) = MatrixFilled_T;
                end
            end
            
        end
        
        function [ ] = KmeansClusterNumSelection(theSignalAnalysis,TrainSet,TestNum,SumplotStates,SilhouetteStates)
            %[] = KmeansClusterNumSelection(TrainngSet,TestNum,SumplotStates,SilhouetteStates)
            %   Find the kmeans clusters numbers by the elbow plots and
            if nargin == 1
                TestNum = 15;
                SumplotStates = 'Y';
                SilhouetteStates = 'Y';
            end
            if nargin == 2
                SumplotStates = 'Y';
                SilhouetteStates = 'Y';
            end
            
            if size(TrainSet,1) < TestNum
                TestNum = size(TrainSet,1) - 1;
            end
            if TestNum == 0
                TestNum = 1;
            end
            
            %   Avg DIs plot
            if ismember(SumplotStates,theSignalAnalysis.BinaryY)
                DisSumInte = zeros(TestNum,1);
                for i = 1:TestNum
                    [~,~,sumd,~] = kmeans(TrainSet,i,'MaxIter',1000,'Replicates',10);
                    % maybe better the avarage distance, but the ideas are similar
                    DisSumInte(i) = sum(sumd)/size(TrainSet,1);
                end
                figure()
                plot([1:TestNum],DisSumInte,'-s','LineWidth',2);
                title('Avg Sum Dis plot');
                ylabel('Avg Sum Dis');
                xlabel('Numbers of clusters');
            end
            %   SilhouetteStates plots
            if ismember(SilhouetteStates,theSignalAnalysis.BinaryY)
                SliplotsNum = ceil(TestNum/4);
                figure;
                for i = 1:4
                    for j = SliplotsNum*(i-1) +1 : SliplotsNum*i
                        idx = kmeans(TrainSet, j,'Replicates',10);
                        subplot(4,SliplotsNum,j)
                        [si,~] = silhouette(gather(TrainSet),gather(idx));
                        title(sprintf('Centroids: %.2g',j))
                        subtitle(sprintf('Mean silhouette value: %g',mean(si)))
                        xlim([-1 1])
                    end
                end
            end
        end
        
        function [X_norm, mu, sigma] = featureNormalize(X)
            %[X_norm, mu, sigma] = featureNormalize(X)
            %FEATURENORMALIZE Normalizes the features in X
            %   FEATURENORMALIZE(X) returns a normalized version of X where
            %   the mean value of each feature is 0 and the standard deviation
            %   is 1.
            
            mu = mean(X);
            X_norm = bsxfun(@minus, X, mu);
            
            sigma = std(X_norm);
            X_norm = bsxfun(@rdivide, X_norm, sigma);
            
        end
        
        function X = NormReconstruct(Xnorm,mu,sigma)
            %X = NormReconstruct(Xnorm, mu,sigma)
            %   Reconstruct the normolized data
            X =  Xnorm * sigma + mu;
        end
        
        function TemplateRawGe = RoughTemplateGe(theSignalAnalysis,TraningSignal,ClusterFeatures,ClusterNum,SignalIntervalIdx)
            TemplateRawGe = cell(ClusterNum,1);
            for i = 1:ClusterNum
                SignInte = TraningSignal(SignalIntervalIdx == i,:);
                %if there is only one signal in this cluster
                if size(SignInte,1) == 1
                    TemplateRawGe{i}  =  nonzeros(SignInte);
                    continue
                else
                    %when it's unmerged interval, the averaging will shift
                    %the less wide peaks to the postion of the widest peak, to
                    %get better templates
                    FeatureInte = ClusterFeatures(SignalIntervalIdx == i,:);
                    [WidestInteLen,WidestIdx] = max(FeatureInte.TimeDuration);
                    %Loc +1 is the peak loc ,6's is the relative loc
                    WidestLoc = round(WidestInteLen*FeatureInte.RelativePeakLoc(WidestIdx)/(1/gather(theSignalAnalysis.Test_freq))) + 1 ;
                    for k = 1:1:size(SignInte,1)
                        if k == WidestIdx
                            continue
                        end
                        PkLoc = round(FeatureInte.TimeDuration(k)*FeatureInte.RelativePeakLoc(k)/(1/theSignalAnalysis.Test_freq)) + 1;
                        SignInte(k,:) = circshift(SignInte(k,:),abs(WidestLoc - PkLoc));
                        
                        % non zeors mean of each interval signal
                        TemplateRawGe{i} = nonzeros(mean(SignInte));
                    end
                end
            end
        end
        %
        function [Regulatedtemplates,RegulatedtemplatesNum] = templatesRegu(theSignalAnalysis,UnRegutemplates,StacountInte,ManSelectedNum)
            %Regulatedtemplates = templatesRegu(theSignalAnalysis,UnRegutemplates,StacountInte,ManSelectedNum)
            %Func to regularize the raw template, take the two side min and
            %one max as the peak templates
            switch isempty(ManSelectedNum)
                case true
                    % find the percentage is over 10%
                    RegulatedtemplatesNum = find(round(StacountInte(:,2)*100) >= 10);
                case false
                    RegulatedtemplatesNum = ManSelectedNum;
            end
            %
            Regulatedtemplates = cell(length(RegulatedtemplatesNum),1);
            for i = 1:length(RegulatedtemplatesNum)
                %                 FeatureInte = theSignalAnalysis.FeatureCentroidReNorm(RegulatedtemplatesNum(i),:);
                Templaterawselected = UnRegutemplates{RegulatedtemplatesNum(i)};
                %exam the peak in the template intervals,a peak is defined
                %as the min in the left interval till max till the right
                %min
                switch theSignalAnalysis.Reactionstates
                    case 'Oxi'
                        [pks,pklocInte,~,~] = findpeaks(Templaterawselected,"WidthReference","halfheight");
                        [maxpksval,maxidx] = max(pks);
                        pklocInte = pklocInte(maxidx);
                        [~,leftMidIdx] = min(abs(maxpksval/2 - Templaterawselected(1:pklocInte)));
                        [~,RightMidIdx] = min(abs(maxpksval/2 - Templaterawselected(pklocInte:end)));
                        RightMidIdx = RightMidIdx + pklocInte - 1;
                        Rsearch = 0.5;
                        try
                            % in the top Rsearch of the peaks find the peak
                            % nearest the midpoint
                            [~,Leftloc] = findpeaks(-Templaterawselected(1:leftMidIdx),"SortStr","descend");
                            [~,maxidx] = min(abs(Leftloc(1:ceil(size(Leftloc,1)*Rsearch)) - leftMidIdx));
                        catch
                            [Leftpks,Leftloc] = max(-Templaterawselected(1:leftMidIdx));
                            [~,maxidx] = max(Leftpks);
                        end
                        if isempty(Leftloc)
                            [Leftpks,Leftloc] = max(-Templaterawselected(1:leftMidIdx));
                            [~,maxidx] = max(Leftpks);
                        end
                        NewLeftidx = Leftloc(maxidx);
                        try
                            [~,Rightloc] = findpeaks(-Templaterawselected(RightMidIdx:end),"SortStr","descend");
                            [~,maxidx] = min(abs(Rightloc(1:ceil(size(Rightloc,1)*Rsearch)) - 1));
                        catch
                            [Rightpks,Rightloc] = max(-Templaterawselected(RightMidIdx:end));
                            [~,maxidx] = max(Rightpks);
                        end
                        if isempty(Rightloc)
                            [Rightpks,Rightloc] = max(-Templaterawselected(RightMidIdx:end));
                            [~,maxidx] = max(Rightpks);
                        end
                        NewRightidx = Rightloc(maxidx) + RightMidIdx - 1;
                    case 'Red'
                        [pks,pklocInte,~,~] = findpeaks(-Templaterawselected,"WidthReference","halfheight");
                        [maxpksval,maxidx] = max(pks);
                        pklocInte = pklocInte(maxidx);
                        [~,leftMidIdx] = min(abs(Templaterawselected(1:pklocInte) - (-maxpksval/2) ));
                        [~,RightMidIdx] = min(abs(Templaterawselected(pklocInte:end) - (-maxpksval/2)));
                        RightMidIdx = RightMidIdx + pklocInte - 1;
                        try
                            % in the top 10% of the peaks find the widest
                            [~,Leftloc] = findpeaks(Templaterawselected(1:leftMidIdx),"SortStr","descend");
                            [~,maxidx] = min(abs(Leftloc(1:ceil(size(Leftloc,1)*Rsearch)) - leftMidIdx));
                        catch
                            [Leftpks,Leftloc] = max(Templaterawselected(1:leftMidIdx));
                            [~,maxidx] = max(Leftpks);
                        end
                        if isempty(Leftloc)
                            [Leftpks,Leftloc] = max(Templaterawselected(1:leftMidIdx));
                            [~,maxidx] = max(Leftpks);
                        end
                        NewLeftidx = Leftloc(maxidx);
                        try
                            [~,Rightloc] = findpeaks(Templaterawselected(RightMidIdx:end),"SortStr","descend");
                            [~,maxidx] = min(abs(Rightloc(1:ceil(size(Rightloc,1)*Rsearch)) - 1));
                        catch
                            [Rightpks,Rightloc] = max(Templaterawselected(RightMidIdx:end));
                            [~,maxidx] = max(Rightpks);
                        end
                        if isempty(Rightloc)
                            [Rightpks,Rightloc] = max(Templaterawselected(RightMidIdx:end));
                            [~,maxidx] = max(Rightpks);
                        end
                        NewRightidx = Rightloc(maxidx) + RightMidIdx - 1;
                end
                Regulatedtemplates{i} = Templaterawselected(NewLeftidx:NewRightidx);
            end
        end
    end
end