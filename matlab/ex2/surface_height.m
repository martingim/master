function [a] = surface_height(run_number)
%%
%run_number = 1
image_names
height_config
close all
folder = day_folder(run_number);
f = frequency(run_number); %The frequency of the run
[omega,T,K,LAMBDA,CP,CG] = wparam(f, 0.33); %calculate the wavenumber and phase velocity
sensor_distance = 0.4; %the dstance between the sensors in meters(should maybe change this to a vector)
%start and stop frame for measuring the elevation of the surface
surface_frame_start = surface_start_stop(run_number, 1);
surface_frame_end = surface_start_stop(run_number, 2);
surface_frame_end = surface_frame_end + rem(surface_frame_end-surface_frame_start+1, 2);%to simplify the fft

sensor_samplerate = 100; %Hz
min_separation = 1/f*sensor_samplerate*0.8; %minimum separation when finding the crests and troughs
%read sensor data from file
sensordata = table2array(readtable(folder + sprintf("/run%d.csv", run_number)));

%remove zeros
sensordata = sensordata(:,5:8);
sensordata(sensordata<350|sensordata>850) = NaN; %replace the outliers with NaN

%convert the data from the sensors to measured surface elevation in meters
measured_surface = [sensor1(sensordata(:,1)), sensor2(sensordata(:,2)) sensor3(sensordata(:,3)) sensor4(sensordata(:,4))];
measured_surface = fillmissing(measured_surface,'linear');
%subtract the mean from the surface to get zero mean
surface_mean = mean(measured_surface, 'omitnan');
measured_surface = measured_surface -surface_mean;

%plot the sensor data from all four sensors with the same waves on top of
%each other
sensor_offset = round(sensor_samplerate*sensor_distance/CP)-2; %sensor offset in samples to plot the same waves on top of each other
figure;
hold on 
plot(measured_surface(:,1))
plot(measured_surface(1*sensor_offset:end,2))
plot(measured_surface(2*sensor_offset:end,3))
plot(measured_surface(3*sensor_offset:end,4))
legend('sensor1', 'sensor2', 'sensor3', 'sensor4')

%Find the crests on troughs from the first sensor and plot with the
%measured data
time = (1:size(measured_surface, 1))/100;
figure;
hold on
surf = measured_surface(:,1);
plot(time, measured_surface(:,1))
LMax = islocalmax(measured_surface(:,1), 'MinSeparation', min_separation);
LMin = islocalmin(measured_surface(:,1), 'MinSeparation', min_separation);
plot(time(LMin), surf(LMin), 'x')
plot(time(LMax), surf(LMax), 'x')

%% Lowpass to smooth the sensor data
figure;
hold on
plot(time, measured_surface(:,1))
lowpass_surface = lowpass(measured_surface, 0.1);
plot(time, lowpass_surface(:,1))
legend('raw', 'lowpass')

%% Without lowpass
%find the crests and troughs of the measured surface to find the mean
%amplitude of the waves.
% surf = measured_surface(surface_frame_start:surface_frame_end,:); %change this to use the sensor_offset
% mean_amplitude = zeros(size(surf, 2),1);
% for i=1:size(surf, 2)
%     LMax = islocalmax(surf(:,i), 'MinSeparation', min_separation);
%     LMin = islocalmin(surf(:,i), 'MinSeparation', min_separation);
%     mean_amplitude(i) = (mean(surf(LMax,i))-mean(surf(LMin, i)))/2;
% end
% 
% a = mean(mean_amplitude);

%% With Lowpass
surf = lowpass_surface(surface_frame_start:surface_frame_end,:); %change this to use the sensor_offset
mean_amplitude = zeros(size(surf, 2),1);
sensor_n = 2;
for i=sensor_n %%Now only uses one sensor
    LMax = islocalmax(surf(:,sensor_n), 'MinSeparation', min_separation);
    LMin = islocalmin(surf(:,sensor_n), 'MinSeparation', min_separation);
    mean_amplitude(i) = (mean(surf(LMax,sensor_n))-mean(surf(LMin, sensor_n)))/2;
end

n_pairs = min([size(surf(LMax,sensor_n), 1) size(surf(LMin, sensor_n), 1)]);
crests = 1;
troughs = 1;
amplitudes = zeros(n_pairs*2-1,1);

for i=1:2:n_pairs*2-1
    amplitudes(i) = 1;
    amplitudes(i+1) = 1;
end
troughs = sum(LMin);
crests = sum(LMin);
number_of_pairs = min(sum(LMin), sum(LMax));

i = 1;
amplitudes = [];
crest_index = 0;
trough_index = 0;
while i<=size(LMin, 1)
    % check if it is a crest or trough
    if LMax(i) == 1
        crest_index = i;
        if trough_index ~= 0
            amplitudes = [amplitudes surf(crest_index, sensor_n)-surf(trough_index, sensor_n)];
        end
    elseif LMin(i) == 1
        trough_index = i;
        if crest_index  ~= 0
            amplitudes = [amplitudes surf(crest_index, sensor_n)-surf(trough_index, sensor_n)];
        end
    end
    %increment the index
    i = i+1;
end


figure;
plot(amplitudes)
title(sprintf('the amplitudes of the waves passing sensor %d', sensor_n));
a = mean(amplitudes)/2;
%msg = sprintf('a=%f, a_lowpass=%f, Diff/a=%f', a, a_low, abs(a-a_low)/a);
%disp(msg);

%% plot the measured surface with theoretical linear and non linear
t = 0:0.001:2;
eta = a*cos(-omega*t); %the linear surface
eta_non_linear = a*cos(-omega*t) +  0.5*a^2*K*cos(2*(-omega*t)) + 3/8*K^2*a^3*cos(3*(-omega*t)); %the non-linear surface

%find the first crest to start plotting from
crests = islocalmax(surf(:,1), 'MinSeparation', 10);
crest_found = -1;
while crest_found<0
    [~, start_idx] = max(crests); %the index to start plotting the measured surface from
    
        
    if surf(start_idx)<0.9*a
        crests(start_idx) = 0;
    else
        crest_found=1;
    end
end

%plot the surfaces
figure
hold on
plot(t, eta)
plot(t, eta_non_linear)
t_measured = 0:0.01:2;
plot(t_measured, surf(start_idx:start_idx+size(t_measured, 2)-1,1)')
legend('linear', 'non-linear', 'measured')


%% FFT of the measured surface
X = fillmissing(surf(:,2),'linear');
Y = fft(X); 

L = size(surf, 1);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = 100;
frequencies = Fs*(0:(L/2))/L;
amp =  max(P1); %alpitude of the first harmonic
figure;
hold on
plot(frequencies, P1, 'LineWidth', 3)
plot(f, 0, 'x')
plot(2*f, 0, 'x')
plot(3*f, 0, 'x')
legend('fft of the measured surface', 'f', '2f', '3f')
xlabel('frequency[Hz]')
xlim([0 f*3.2])

%% save the results to params.mat
load params.mat params
p = params(run_number);
p('amplitude_first_harmonic') = amp;
p('a') = a;
p('std_a') = std(amplitudes);
%p('a_low') = a_low;
params(run_number) = p;
save('params.mat', 'params')
end