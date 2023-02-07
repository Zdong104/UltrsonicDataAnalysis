clc
clear all
close all
% process data
sample_cond = {};
for cnt = 1:32
    if cnt < 10
        FileName = strcat('Sample_Condensed_5MHz/Sample_Condensed_5MHz_0',num2str(cnt),".mat");
    else
        FileName = strcat('Sample_Condensed_5MHz/Sample_Condensed_5MHz_', num2str(cnt),".mat");
    end
    sample_cond{cnt} = load(FileName);
    A = [];
    A = sample_cond{cnt}.A;
    for interval = 1:size(A)
        if A(interval)> 4.5*10^(-3)
            x1(cnt) = interval*sample_cond{cnt}.Tinterval;
            break
        end
    end
    dt_w(cnt) = sample_cond{cnt}.Tinterval;
    TOF_r(cnt) = dt_w(cnt) * x1(cnt);
end

Water_cond = {};
for cnt = 1:32
    if cnt < 10
        FileName = strcat('Water_Condensed_5MHz/Water_Condensed_5MHz_0',num2str(cnt),".mat");
    else
        FileName = strcat('Water_Condensed_5MHz/Water_Condensed_5MHz_', num2str(cnt),".mat");
    end
    Water_cond{cnt} = load(FileName);
    A = [];
    A = Water_cond{cnt}.A;
    for interval = 1:size(A)
        if A(interval)> 4.5*10^(-3)
            x1(cnt) = interval*Water_cond{cnt}.Tinterval;
            break
        end
    end
    dt_s(cnt) = Water_cond{cnt}.Tinterval;
    TOF_s(cnt) = dt_s(cnt) * x1(cnt);
end

%%Calculate the speed of material
L= 0.023; % meter
c_r= 1500; % speed of sound in water meter per second

for count = 1:32 % for each data of each trial
    c_s(count)= L/(TOF_s(count) - TOF_r(count)+L/c_r);
end


% additional property
mass_condesned = 0.0046; % kilogram
volume_condensed = 0.023*0.023*0.024; % meter unit
density_condensed = mass_condesned/volume_condensed;
density_water = 1000; % kg/m^3
for count = 1:32
Z_s(count)= density_condensed * c_s(count);
Z_r= density_water*c_r;
T(count)=(2*Z_s(count))/(Z_s(count)+Z_r);
end

% % frequency-dependentattenuation

% a_sub = 20*log10(exp(1))*(1/L(log(A_r/A_s)-log(T)))

%% save file 
%writematrix(TOF_r,'Trial time of flight.txt')

save('Tiral2.mat','TOF_r')