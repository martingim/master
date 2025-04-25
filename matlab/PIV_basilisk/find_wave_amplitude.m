% %choose which day to use the results from 
%%
lab_folder = "Lab_results/24_04_12/";
lab_surface_probes = ["f1425_a0_308.csv";
            "f1425_a0_308_r2.csv";
            "f1425_a0_308_r3.csv";
            "f1425_a0_616_r4.csv";  
            "f1425_a0_154_r5.csv"];


%%
lab_folder = "Lab_results/24_09_18/";
lab_surface_probes = ["f1425_a308_r1.csv";
            "f1425_a308_r2.csv";
            "f1425_a308_r3.csv";
            "f1425_a616_r4.csv";
            "f1425_a616_r5.csv";
            "f1425_a616_r6.csv";
            "f1425_a154_r7.csv";
            "f1425_a154_r8.csv";
            "f1425_a154_r9.csv"];

%%
close all
sensor_samplerate = 100; %Hz
f = 1.425;
min_separation = 1/f*sensor_samplerate*0.8; %minimum separation when finding the crests and troughs
start_idx = 3900;
end_idx = 7000;


run_number = 1;
for run_number=1:5
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
        plot(t(LMax(:,i)), t(LMax(:,i))*0+mean(wave_amplitudes{i})-std(wave_amplitudes{i}), '--', Color=colors(i))
        plot(t(LMax(:,i)), t(LMax(:,i))*0+mean(wave_amplitudes{i})+std(wave_amplitudes{i}), '--', Color=colors(i))
    end
    legend();
    title(sprintf("run number %d", run_number))
    std(wave_amplitudes{4})
    % 
    % plot(t(LMax), smoothed_data(LMax, 1), 'x');
end
