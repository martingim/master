% %choose which day to use the results from 
%%
lab_folder = "Lab_results/24_04_12/";
lab_surface_probes = ["f1425_a0_308.csv";
            "f1425_a0_308_r2.csv";
            "f1425_a0_308_r3.csv";
            "f1425_a0_616_r4.csv";  
            "f1425_a0_154_r5.csv"];
sensor_locations = [8 10.05 10.75 11.50];
start_indices = [3900 3900 3900 3900 3900];
end_indices = [7000 7000 7000 7000 7000];
%%
lab_folder = "Lab_results/24_09_18/";
lab_surface_probes = ["f1425_a308_r1.csv";
            "f1425_a308_r2.csv";
            "f1425_a308_r3.csv";
            %"f1425_a616_r4.csv";
            %"f1425_a616_r5.csv";
            "f1425_a616_r6.csv";
            "f1425_a154_r7.csv";
            "f1425_a154_r8.csv";
            "f1425_a154_r9.csv"];
sensor_locations = [1.5 10.05 10.75 11.50];
start_indices = [3900 3900 3900 3900 3900];
end_indices = [7000 7000 7000 7000 7000];

%%
close all
sensor_samplerate = 100; %Hz
f = 1.425;
min_separation = 1/f*sensor_samplerate*0.8; %minimum separation when finding the crests and troughs

amplitudes = {};

run_number = 1;
for run_number=1:5
    start_idx = start_indices(run_number);
    end_idx = end_indices(run_number);

    sensor_data = load(append(lab_folder, lab_surface_probes(run_number)));
    sensor_data = -sensor_data(:,3:6);
    smoothed_data = lowpass(sensor_data, .7);
    t = (1:size(sensor_data, 1))/sensor_samplerate;
    sensor_data = sensor_data(start_idx:end_idx,:);
    smoothed_data = smoothed_data(start_idx:end_idx,:);
    t = t(start_idx:end_idx);
    figure
    hold on
    for i=1:4
        plot(t, sensor_data(:,i), '--', 'DisplayName',sprintf("sensor %d raw data", i))
        plot(t, smoothed_data(:,i), 'DisplayName',sprintf("sensor %d smoothed", i))
    end
    legend();
    LMax = islocalmax(smoothed_data, 1, 'MinSeparation',min_separation);
    LMin = islocalmin(smoothed_data, 1, 'MinSeparation',min_separation);
    for i=1:4
        n_maxima = size(smoothed_data(LMax(:,i),i),1);
        n_minima = size(smoothed_data(LMin(:,i),i),1);
        if n_maxima>n_minima
            LMax(find(LMax(:,i),1, 'last'),i) = 0;
        elseif n_maxima<n_minima
            LMin(find(LMin(:,i),1, 'last'),i) = 0;
        end
    end
    
    wave_amplitudes = {};
    
    figure
    hold on
    for i=1:4
        run_number
        size(smoothed_data(LMax(:,i),i))
        size(smoothed_data(LMin(:,i),i))
        wave_amplitudes{i} = (smoothed_data(LMax(:,i),i)-smoothed_data(LMin(:,i),i))/2;
        plot(t, smoothed_data(:,i), 'DisplayName',sprintf("sensor %d", i))
        plot(t(LMax(:,i)), smoothed_data(LMax(:,i),i), "x")
        plot(t(LMin(:,i)), smoothed_data(LMin(:,i),i), "x")
        % mean_amplitude(i) = (mean(smoothed_data(LMax,i))-mean(smoothed_data(LMin, i)))/2;
    end
    legend();
    figure
    hold on
    colors = ["red" "blue" "green" "black"];
    for i=1:4
        plot(t(LMax(:,i)),wave_amplitudes{i},'x','DisplayName',sprintf("sensor %d", i), Color=colors(i))
        plot(t(LMax(:,i)), t(LMax(:,i))*0+mean(wave_amplitudes{i}), 'DisplayName', sprintf("mean amplitude sensor %d", i),Color=colors(i))
        plot(t(LMax(:,i)), t(LMax(:,i))*0+mean(wave_amplitudes{i})-std(wave_amplitudes{i}), '--', 'DisplayName', sprintf("sensor %d amplitude + %s", i, '\sigma'), Color=colors(i))
        plot(t(LMax(:,i)), t(LMax(:,i))*0+mean(wave_amplitudes{i})+std(wave_amplitudes{i}), '--', 'DisplayName', sprintf('sensor %d amplitude - %s', i, '\sigma'), Color=colors(i))
    end
    legend();
    title(sprintf("run number %d", run_number))


    run_amplitudes = [];    
    for i=1:4
        run_amplitudes = [run_amplitudes mean(wave_amplitudes{i})];
    end
    amplitudes{run_number} = run_amplitudes;
end

figure
hold on
for run_number =5
    plot(sensor_locations, amplitudes{run_number}/amplitudes{run_number}(1), '-x','DisplayName',sprintf("run number %d", run_number))
end
xlabel('sensor distance from piston [m]')
ylabel("measured amplitude")
title('measured amplitudes different sensors')
legend()

%% run number 1
lab_run_number = 1;
basilisk_folders = [];
titles = [];
runtimes = [];
start_indices =[];
end_indices =[];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL11_0/"]; titles = [titles "LEVEL 11"]; runtimes=[runtimes 843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL11_1/"]; titles = [titles "LEVEL 11+1"];runtimes=[runtimes 1792];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL11_2/"]; titles = [titles "LEVEL 11+2"];runtimes=[runtimes 4308];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL12_0/"]; titles = [titles "LEVEL 12+0"];runtimes=[runtimes 4089];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL12_1/"]; titles = [titles "LEVEL 12+1"];runtimes=[runtimes 4089];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL13_0/"]; titles = [titles "LEVEL 13+0"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL11/"]; titles = [titles "boundary LEVEL 11"];runtimes=[runtimes 350];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL12/"]; titles = [titles "boundary LEVEL 12"];runtimes=[runtimes 6843];
%basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL13/"]; titles = [titles "boundary LEVEL 13"];runtimes=[runtimes 6843];
start_indices = [800 2400 2400 2400];
end_indices = [2500 2500 2500 2500 ];


%plot the basilisk results
Fs = 100;
f = 1.425;
min_separation = 1/f*Fs*0.8; %minimum separation when finding the crests and troughs

basilisk_wave_amplitudes = {};
for i=1:size(basilisk_folders, 2)
    
    basilisk_folder = basilisk_folders(i);
    %plot the basilisk results
    surface_probes = readtable(append(basilisk_folder, "surface_probes.csv"), 'ReadVariableNames',true, 'VariableNamingRule','preserve');
    probe_names = surface_probes.Properties.VariableNames;
    probe_locations = str2double(surface_probes.Properties.VariableNames);
    surface_probes = table2array(surface_probes);

    basilisk_run_amplitudes = [];
    figure;
    hold on
    for sensor=1:4
        basilisk_probe_number = sensor+1;
        basilisk_surface = surface_probes(start_indices(sensor):end_indices(sensor),basilisk_probe_number);
        t = 0:1/Fs:(size(surface_probes(:,basilisk_probe_number), 1)-1)/Fs;
        t = t(start_indices(sensor):end_indices(sensor));
        LMax = islocalmax(basilisk_surface, 1, 'MinSeparation',min_separation);
        LMin = islocalmin(basilisk_surface, 1, 'MinSeparation',min_separation);
        n_maxima = size(basilisk_surface(LMax),1);
        n_minima = size(basilisk_surface(LMin),1);
        if n_maxima>n_minima
            LMax(find(LMax,1, 'last')) = 0;
        elseif n_maxima<n_minima
            LMin(find(LMin,1, 'last')) = 0;
        end
        basilisk_run_amplitude = (basilisk_surface(LMax)-basilisk_surface(LMin))/2;
        basilisk_run_amplitudes = [basilisk_run_amplitudes mean(basilisk_run_amplitude)];
        plot(t(LMax), basilisk_run_amplitude,'x', "DisplayName",sprintf("probe %d", sensor))
    end
    legend()
    basilisk_wave_amplitudes{i} = basilisk_run_amplitudes;
end

%plot the evolution of the amplitude
close all
figure
hold on
for i=1:size(basilisk_folders,2)
    plot(probe_locations(2:5), basilisk_wave_amplitudes{i}, 'DisplayName',titles(i))
end

plot(sensor_locations, amplitudes{1}, '-x','DisplayName',sprintf("Lab experiment run number %d", run_number))
legend()
