clc
clear all
close all

% sampleNum = 4;
sampleThickness = 0.053;
cuttingTime = 3; %mu s

%%
% temperatures = [1, 2, 3, 4, 5, 8, 9, 10];
% 
% for tempCnt = 1 : length(temperatures)
%     
%    for expCnt = 1:3
       close all;
       %% Medium Signal
        expDirectorySponge = strcat('9_15/Wax/Temp',num2str(8));
        FolderNameMedium = strcat('20170921-000',num2str(3));

        for cnt = 1 : 32

            if cnt < 10
                FileName = strcat(expDirectorySponge,'/',FolderNameMedium,'/',FolderNameMedium,'_0',num2str(cnt),'.csv');
            else
                FileName = strcat(expDirectorySponge,'/',FolderNameMedium,'/',FolderNameMedium,'_',num2str(cnt),'.csv');
            end

            [text1,text2] = textread(FileName,'%s %s','delimiter',',');
            voltageIndicatorMedium = text2(2,:);
            text1 = (text1(4:end));
            text2 = (text2(4:end));

            indexC = strfind(text2,'Infinity');
            index = find(not(cellfun('isempty', indexC)));

            if ~isempty(index)
                text2(index) = {num2str(0)};
            end

            xDataS = sprintf('%s*',text1{:});
            xData = sscanf(xDataS,'%f*');

            yDataS = sprintf('%s*',text2{:});
            yData = sscanf(yDataS,'%f*');    

        %     signal_Sponge(:,cnt*2-1:cnt*2) = csvread(FileName,2,0);
            signal_Sponge(:,cnt*2-1:cnt*2) = [xData,yData];
            if strcmp(voltageIndicatorMedium,'(mV)')
                ampSponge(:,cnt) = signal_Sponge(:,cnt*2)*10^(-3);
            else
                ampSponge(:,cnt) = signal_Sponge(:,cnt*2);
            end


        end

        ampAvg_Medium(1,:) = signal_Sponge(:,1);
        ampAvg_Medium(2,:) = mean(transpose(ampSponge));
        % I = find (signal_Sponge(:,1) >= 0);temp(1,:) = signal_Sponge(I,1);temp(2,:) = ampAvg_Sponge(I);
        % I = find (temp(1,:) <= 5);amgAvgCutSponge(1,:) = temp(1,I);amgAvgCutSponge(2,:) = temp(2,I);
        % ampAvg_Sponge = amgAvgCutSponge;

        %% Water Signal

        expDirectoryWater =  expDirectorySponge;
        FolderNameWater = strcat('Water');


        for cnt = 1 : 27

            if cnt < 10
                FileName = strcat(expDirectoryWater,'/',FolderNameWater,'/',FolderNameWater,'_0',num2str(cnt),'.csv');
            else
                FileName = strcat(expDirectorySponge,'/',FolderNameWater,'/',FolderNameWater,'_',num2str(cnt),'.csv');
            end

            [text1,text2] = textread(FileName,'%s %s','delimiter',',');
            voltageIndicatorWater = text2(2,:);

            text1 = (text1(4:end));
            text2 = (text2(4:end));

            indexC = strfind(text2,'Infinity');
            index = find(not(cellfun('isempty', indexC)));

            if ~isempty(index)
                text2(index) = {num2str(0)};
            end

            xDataS = sprintf('%s*',text1{:});
            xData = sscanf(xDataS,'%f*');

            yDataS = sprintf('%s*',text2{:});
            yData = sscanf(yDataS,'%f*');    

            signal_Water(:,cnt*2-1:cnt*2) = [xData,yData];
             if strcmp(voltageIndicatorWater,'(mV)')
                ampWater(:,cnt) = signal_Water(:,cnt*2)*10^(-3);
            else
                ampWater(:,cnt) = signal_Water(:,cnt*2);
            end


%             ampWater(:,cnt) = signal_Water(:,cnt*2);

        end
        ampAvg_Water(1,:) = signal_Water(:,1);
        ampAvg_Water(2,:) = mean(transpose(ampWater));
        % I = find (signal_Water(:,1) >= 0);temp(1,:) = signal_Water(I,1);temp(2,:) = ampAvg_Water(I);
        % I = find (temp(1,:) <= 5);amgAvgCutWater(1,:) = temp(1,I);amgAvgCutWater(2,:) = temp(2,I);
        % ampAvg_Water = amgAvgCutWater;

        %% Cut the signal

        cuttingIndexMedium = find ( ampAvg_Water(1,:) < cuttingTime & ampAvg_Water(1,:) > 0);
        cutSignalMedium(1,:) = ampAvg_Medium(1,cuttingIndexMedium);
        cutSignalMedium(2,:) = ampAvg_Medium(2,cuttingIndexMedium);

        cuttingIndexWater = find ( ampAvg_Water(1,:) < cuttingTime & ampAvg_Water(1,:) > 0);
        cutSignalWater(1,:) = ampAvg_Water(1,cuttingIndexWater);
        cutSignalWater(2,:) = ampAvg_Water(2,cuttingIndexWater);

        % No Cut
        % cutSignalWater(1,:) = ampAvg_Water(1,:);
        % cutSignalWater(2,:) = ampAvg_Water(2,:);

        %% FFT 
        addpath('/home/Omid/Projects/omidMatlabToolBox'); %ToolBox for FFT Calculating

        lengthCoeff = 1;
        samplingFrequency = 1 /  ( cutSignalMedium(1,2) - cutSignalMedium(1,1) );

        % fftLength_Bone = length(ampAvg_Sponge(2,:)) * lengthCoeff;
        % FFT_Bone = fft(ampAvg_Sponge(2,:),fftLength_Bone);
        % FFT_Bone = FFT_Bone(1:fftLength_Bone/2);
        % magFFT_Bone = abs(FFT_Bone);
        % f_Bone = (0:fftLength_Bone/2-1) * samplingFrequency / fftLength_Bone;

        FFT_Bone = OmidFFTCaclculator(cutSignalMedium(2,:),samplingFrequency,5);
        magFFT_Bone = FFT_Bone.magnitude;
        f_Bone = FFT_Bone.frequency;

        FFT_Water = OmidFFTCaclculator(cutSignalWater(2,:),samplingFrequency,5);
        magFFT_Water = FFT_Water.magnitude;
        f_Water = FFT_Water.frequency;

        % fftLength_Water = length(cutSignalWater(2,:)) * lengthCoeff;
        % FFT_Water = fft(ampAvg_Water(2,:),fftLength_Water);
        % FFT_Water = FFT_Water(1:fftLength_Water/2);
        % magFFT_Water = abs(FFT_Water);
        % f_Water = (0:fftLength_Water/2-1) * samplingFrequency / fftLength_Water;

        %% Curve fitting
        windowLength_Bone = 500;
        g_Bone = gausswin (windowLength_Bone);
        g_Bone = g_Bone / sum(g_Bone);


        windowLength_Water = 500;
        g_Water = gausswin (windowLength_Water);
        g_Water = g_Water / sum(g_Water);

        magFFT_Bone_Smooth = conv(magFFT_Bone, g_Bone, 'same');
        magFFT_Water_Smooth = conv(magFFT_Water, g_Water, 'same');

        %% Plotting

        fig1 = figure('Name','Medium Signal');
        subplot(2,1,1)
        plot(cutSignalMedium(1,:),cutSignalMedium(2,:))
        xlabel('Time(\mu s)'); ylabel('Amplitude');
        subplot(2,1,2)
        plot(f_Bone,magFFT_Bone)
        xlabel('Frequency(kHz)'); ylabel('Amplitude');
        hold on
        plot(f_Bone,magFFT_Bone_Smooth,'r')
        hold off
        % xlim([0 5])

        fig2 = figure('Name','Water Signal');
        subplot(2,1,1)
        plot(cutSignalWater(1,:),cutSignalWater(2,:))
        subplot(2,1,2)
        plot(f_Water,magFFT_Water)
        hold on
        plot(f_Water,magFFT_Water_Smooth,'r')
        hold off
        % xlim([0 5])

        fig3 = figure('Name','Comparison');
        subplot(2,1,1)
        plot(f_Water,magFFT_Water_Smooth)
        hold on
        plot(f_Bone,magFFT_Bone_Smooth,'r')
        % xlim([0 5])
        subplot(2,1,2)
        plot(f_Water,magFFT_Water_Smooth/max(magFFT_Water_Smooth))
        hold on
        plot(f_Bone,magFFT_Bone_Smooth/max(magFFT_Bone_Smooth),'r')
        % xlim([0 5])

        %% Attenuation
        % alpha = - 20 * log10(exp(1)) * (1/sampleThickness) * log(magFFT_Bone_Smooth./magFFT_Water_Smooth);
        alpha = - (1/sampleThickness) * log(magFFT_Bone_Smooth./magFFT_Water_Smooth);

        % Find the attenuation value at central frequency 1 MHz

        I = find( (870 <= f_Bone) & (f_Bone<= 1280) );f_Bone = f_Bone(I);alpha=alpha(I);

        fc = 1000;
        I_cent = find( abs(f_Bone-fc) == min( abs(f_Bone-fc) ) );
        % centralFreq = f_Bone(I_cent)
        % centralAtten = alpha(I_cent)

        polyfitCoeff = polyfit(f_Bone,alpha,1);
        alpha_Lin = polyval(polyfitCoeff,f_Bone);
        centralAtten = alpha_Lin(I_cent);
        BUA = polyfitCoeff(1);

        fig4 = figure('Name','Attenuation');
        plot(f_Bone,alpha); hold on; plot(f_Bone,alpha_Lin,'r')
        xlabel('Frequency (kHz)')
        ylabel('Attenuation (1/cm)')
        % ylim([0 50])
        % title('Sample 1 - Experiment 3')
        saveas(fig1,strcat(expDirectorySponge,'/',FolderNameMedium,'/MediumSignal.fig'))
        saveas(fig2,strcat(expDirectorySponge,'/',FolderNameMedium,'/WaterSignal.fig'))
        saveas(fig3,strcat(expDirectorySponge,'/',FolderNameMedium,'/Comparison.fig'))
        saveas(fig4,strcat(expDirectorySponge,'/',FolderNameMedium,'/Attenuation.fig'))
        %% Speed of Sound

        % I_segment1 = find( (0<=cutSignalMedium(1,:)) & (cutSignalMedium(1,:)<=0.02) ); segment1 = cutSignalMedium(:,I_segment1);
        % I_segment2 = find( (0.02<=cutSignalMedium(1,:)) & (cutSignalMedium(1,:)<=0.04) ); segment2 = cutSignalMedium(:,I_segment2);

        threshold = max(abs(cutSignalMedium(2,:))) / 5.3;

        % I_TOF = find ( min(abs(abs(segment2(2,:))-threshold)) == abs(abs(segment2(2,:))-threshold) );

        segment = abs(cutSignalMedium(2,:)) - threshold;
        I_TOF = find(segment>0);
        TOF = cutSignalMedium(1,I_TOF(1));



        % I_max = find(max(cutSignalMedium(2,:)) == cutSignalMedium(2,:));
        % threshold = abs(cutSignalMedium(2,I_max)) / 14.3;
        % I_TOF = find ( min(abs(abs(cutSignalMedium(2,:))-threshold)) == abs(abs(cutSignalMedium(2,:))-threshold) );

        % TOF = cutSignalMedium(1,I_TOF); %mu s
        SOS = sampleThickness/ TOF;
        save(strcat(expDirectorySponge,'/',FolderNameMedium,'/Parameters.mat'),'centralAtten','BUA','SOS','TOF')

%        
%    end
%     
% end


