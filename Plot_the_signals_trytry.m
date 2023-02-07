clc
clear all
close all
%% Load Sample Now
FileName = load(['2.25 condensed.mat']);
A = [];
    A = FileName.A;
    for interval_Sample = 1:size(A)
        if A(interval_Sample)> 0.002
            Signal_Sample = interval_Sample;
            break
        end
    end
figure(1)
plot(A)% meaningless plot since it does not have any units
title("Raw data Collected from Probe for 5MHz Loose Sample")
saveas(figure(1),'Raw Data for 5MHz Loose Sample.jpg')
%%
interval_end_Sample = interval_Sample+100
interval_start_Sample = interval_Sample-20
Signal_Sample = A(interval_start_Sample:interval_end_Sample,:); %water only received Signal_Samples
%Signal_Sample = A(156200:157000,:); % received Signal_Sample with sample in the Loose
figure(2)
plot(Signal_Sample) %meaningless plot!!
dt= FileName.Tinterval;% in milli seconds since the units chosen on picoscope was ms!
fs=1/dt;% it is in kHz
desired_freq = 1*2250; % This is the desired frequency in kHz which is 1 MHz
figure(3)
time_vec_s = 0:dt*10^6:length(Signal_Sample)*dt*10^6;
plot(time_vec_s(2:end)',Signal_Sample)
xlabel('Time Domain, \mu s')
ylabel('Recived Signal, Volts')

figure(4) % find TOF_s
time_vec_A= 0:dt*10^6:length(A)*dt*10^6;
plot(time_vec_A(2:end)',A)
xlabel('Time , \mus')
ylabel('Volts')

% Now, I have a plot for my signal in terms of time, taking into account my
% time step 

% here signal is missing, so it will cause problem for finding the fft of
% the signal
% we do the interpolate to find the missing values
% for i=1:length(Signal)
%     Signal(isinf(Signal))=NaN;
% end
% [F,TF] = fillmissing(Signal,'spline');
% plot(F)
%% Now, for the Fourier transform of this:
% Signal=F;
fftSignal = abs(fft(Signal_Sample,10000));
% figure
% plot(fftSignal)

%fs = 1/dt;
freq_vec = 0:fs/length(fftSignal):(fs-fs/length(fftSignal));
figure(5)
plot(freq_vec,fftSignal)
title('Frequency Domained Amplitude for Sample')
xlim([0,2000*desired_freq])
xlabel('Frequency Domain, Hz')
ylabel('Amplitude')
 %% Curve fitting
        windowLength_Sample = interval_end_Sample-interval_start_Sample;
        g_Sample = gausswin (windowLength_Sample);
        g_Sample = g_Sample / sum(g_Sample);

        magFFT_Sample_Smooth = conv(fftSignal, g_Sample, 'same');
       
figure(6)  
plot(freq_vec, magFFT_Sample_Smooth)
title('Curve fitting for Sample')
saveas(figure(6),'Curve fitting for Sample.jpg')
Max_amplitude_Sample = max(magFFT_Sample_Smooth);% this is the max point for amplitude 


%% Load water data now
FileName = load('2.25 water.mat');
A = [];
    A = FileName.A;
    for interval_water = 1:size(A)
        if A(interval_water)> 0.005
            Signal_Water = interval_water;
            break
        end
    end
figure(7)
plot(A)
title('Raw data Collected from Probe for 5MHz Loose Water')
saveas(figure(7),'Raw Data for 5MHz Loose Water.jpg')% meaningless plot since it does not have any units
%%
interval_end_water = interval_water+100
interval_start_water = interval_water-20
Signal_Water = A(interval_start_water:interval_end_water,:); %water only received Signal_Waters
%Signal_Water = A(156200:157000,:); % received Signal_Water with sample in the Loose
figure(8)
plot(Signal_Water) %meaningless plot!!
dt=FileName.Tinterval;% in milli seconds since the units chosen on picoscope was ms!
fs=1/dt;% it is in kHz
desired_freq = 1*2250; % This is the desired frequency in kHz which is 1 MHz
figure(9)
time_vec_w = 0:dt*10^6:length(Signal_Water)*dt*10^6;
plot(time_vec_w(2:end)',Signal_Water)
xlabel('Time Domain, \mu s')
ylabel('Recived Signal, Volts')

figure(10) % find TOF_s
time_vec_A= 0:dt*10^6:length(A)*dt*10^6;
plot(time_vec_A(2:end)',A)
xlabel('Time , \mus')
ylabel('Volts')

% Now, I have a plot for my signal in terms of time, taking into account my
% time step 

% here signal is missing, so it will cause problem for finding the fft of
% the signal
% we do the interpolate to find the missing values
% for i=1:length(Signal)
%     Signal(isinf(Signal))=NaN;
% end
% [F,TF] = fillmissing(Signal,'spline');
% plot(F)
%% Now, for the Fourier transform of this:
% Signal=F;
fftSignal = abs(fft(Signal_Water,10000));
% figure
% plot(fftSignal)

%fs = 1/dt;
freq_vec = 0:fs/length(fftSignal):(fs-fs/length(fftSignal));
figure(11)
plot(freq_vec,fftSignal)
title('Frequency Domained Amplitude for reference')
xlim([0,2000*desired_freq])
xlabel('Frequency Domain, Hz')
ylabel('Amplitude')

 %% Curve fitting
        windowLength_Water =interval_end_water-interval_start_water;
        g_Water = gausswin (windowLength_Water);
        g_Water = g_Water / sum(g_Water);

        magFFT_Water_Smooth = conv(fftSignal, g_Water, 'same');
       
figure(12)  
plot(freq_vec, magFFT_Water_Smooth)
title('Curve fitting for reference')
saveas(figure(12),'Curve fitting for Water.jpg')
Max_amplitude_Water = max(magFFT_Water_Smooth);
% Now, I want to plot the dB of this
% dBfftSignal = 20*log10(fftSignal);
% figure
% plot(freq_vec,dBfftSignal)
% xlim([0,2*desired_freq])
% xlabel('Frequency Domain, MHz')
% ylabel('dB')
% pr=1;
% ps=;
% cs=;
% cr=;
A_s = Max_amplitude_Sample;
A_r = Max_amplitude_Water;% for Max_amplitude for water

density_water = 0.001; % g/mm^3
mass_condesned = 2; % kilogram
L = 25; % mm
c_r= 1500; %mm/mu s
volume_Loose = 26*24*25; % mm unit
density_Loose = mass_condesned/volume_Loose;
TOFs = interval_start_Sample;
TOFr = interval_start_water;
c_s= L/(TOFs-TOFr+(L/c_r)); % speed of sound in water meter per second
Z_r= density_water*c_r ;
Z_s = density_Loose * c_s;
T =(2*Z_s)/(Z_s+Z_r);
alpha_sub=20*log10(exp(1))*((1/L)*(log(A_r/A_s)-log(T)))



