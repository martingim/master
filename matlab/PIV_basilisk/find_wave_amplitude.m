% %choose which day to use the results from 
%%
lab_folder = "Lab_results/24_04_12/";
ex_day = "12/04/24";
run_numbers = 1:3;
lab_surface_probes = ["f1425_a0_308.csv";
            "f1425_a0_308_r2.csv";
            "f1425_a0_308_r3.csv";
            "f1425_a0_616_r4.csv";  
            "f1425_a0_154_r5.csv"];
sensor_locations = [8 10.05 10.75 11.50];
start_indices = [3900 3900 3900 3900 3900;
                3900 3900 3900 3900 3900;
                3900 3900 3900 3900 3900;
                3900 3900 3900 3900 3900;
                3900 3900 3900 3900 3900];
end_indices = [5000 5000 5000 5000 5000;
                5000 5000 5000 5000 5000;
                5000 5000 5000 5000 5000;
                5000 5000 5000 5000 5000;
                5000 5000 5000 5000 5000;];
%%
lab_folder = "Lab_results/24_09_18/";
ex_day = "18/09/24";
run_numbers = [1 2 3 6 7 8 9];
% run_numbers = [7 8 9];
lab_surface_probes = ["f1425_a308_r1.csv";
            "f1425_a308_r2.csv";
            "f1425_a308_r3.csv";
            "f1425_a616_r4.csv";
            "f1425_a616_r5.csv";
            "f1425_a616_r6.csv";
            "f1425_a154_r7.csv";
            "f1425_a154_r8.csv";
            "f1425_a154_r9.csv"];
sensor_locations = [1.5 10.05 10.75 11.50];
start_indices =[820 3000 3000 3000;
                820 3000 3000 3000;
                820 3000 3000 3000;
                800 2600 2800 3000;
                800 2600 2800 3000;
                800 2600 2800 3000;
                800 2600 2800 3000;
                800 2600 2800 3000;
                800 2600 2800 3000];
end_indices = [ 5000 5000 5000 5000;
                5000 5000 5000 5000;
                5000 5000 5000 5000;
                5000 5000 5000 5000;
                5000 5000 5000 5000;
                5000 5000 5000 5000;
                5000 5000 5000 5000;
                5000 5000 5000 5000;
                5000 5000 5000 5000];

%%
% close all
sensor_samplerate = 100; %Hz
f = 1.425;
min_separation = 1/f*sensor_samplerate*0.8; %minimum separation when finding the crests and troughs
plot_results = true;
amplitudes = {};

for run_number=run_numbers
    start_idx = start_indices(run_number,:);
    end_idx = end_indices(run_number,:);

    sensor_data = load(append(lab_folder, lab_surface_probes(run_number)));
    sensor_data = -sensor_data(:,3:6);
    smoothed_data = lowpass(sensor_data, .7);
    t = (1:size(sensor_data, 1))/sensor_samplerate;
    if plot_results
        figure
        hold on
    end
    LMax = {};
    LMin = {};
    surface_elevation  = {};
    smoothed_elevation = {};
    t_sensor = {};
    for sensor=1:4
        surface_elevation{sensor} = sensor_data(start_idx(sensor):end_idx(sensor),sensor);
        smoothed_elevation{sensor} = smoothed_data(start_idx(sensor):end_idx(sensor),sensor);
        t_sensor{sensor} = t(start_idx(sensor):end_idx(sensor));
        if plot_results
            plot(t_sensor{sensor}, surface_elevation{sensor}, '--', 'DisplayName',sprintf("sensor %d raw data", sensor))
            plot(t_sensor{sensor}, smoothed_elevation{sensor}, 'DisplayName',sprintf("sensor %d smoothed", sensor))
        end
        LMax{sensor} = islocalmax(smoothed_elevation{sensor}, 1, 'MinSeparation',min_separation);
        LMin{sensor} = islocalmin(smoothed_elevation{sensor}, 1, 'MinSeparation',min_separation);
        n_maxima = size(smoothed_elevation{sensor}(LMax{sensor}),1);
        n_minima = size(smoothed_elevation{sensor}(LMin{sensor}),1);
        if n_maxima>n_minima
            LMax{sensor}(find(LMax{sensor},1, 'last')) = 0;
        elseif n_maxima<n_minima
            LMin{sensor}(find(LMin{sensor},1, 'last')) = 0;
        end
    end
    if plot_results
        legend();
        figure
        hold on
    end
    
    wave_amplitudes = {};
    for i=1:4
        wave_amplitudes{i} = (smoothed_elevation{i}(LMax{i})-smoothed_elevation{i}(LMin{i}))/2;
        if plot_results
            plot(t_sensor{i}, smoothed_elevation{i}, 'DisplayName',sprintf("sensor %d", i))
            plot(t_sensor{i}(LMax{i}), smoothed_elevation{i}(LMax{i}), "x")
            plot(t_sensor{i}(LMin{i}), smoothed_elevation{i}(LMin{i}), "x")
            % mean_amplitude(i) = (mean(smoothed_data(LMax,i))-mean(smoothed_data(LMin, i)))/2;
        end
    end
    
    if plot_results
        legend();
        figure
        hold on
        colors = ["red" "blue" "green" "black"];
        for i=1:4
            plot(t_sensor{i}(LMax{i}), wave_amplitudes{i},'x','DisplayName',sprintf("sensor %d", i), Color=colors(i))
            plot(t_sensor{i}(LMax{i}), t_sensor{i}(LMax{i})*0+mean(wave_amplitudes{i}), 'DisplayName', sprintf("mean amplitude at x=%.2f", sensor_locations(i)),Color=colors(i))
            plot(t_sensor{i}(LMax{i}), t_sensor{i}(LMax{i})*0+mean(wave_amplitudes{i})-std(wave_amplitudes{i}), '--', 'DisplayName', sprintf("sensor %d amplitude + %s", i, '\sigma'), Color=colors(i))
            plot(t_sensor{i}(LMax{i}), t_sensor{i}(LMax{i})*0+mean(wave_amplitudes{i})+std(wave_amplitudes{i}), '--', 'DisplayName', sprintf('sensor %d amplitude - %s', i, '\sigma'), Color=colors(i))
        end
        legend();
        title(sprintf("run number %d", run_number))
    end

    run_amplitudes = [];    
    for i=1:4
        run_amplitudes = [run_amplitudes mean(wave_amplitudes{i})];
    end
    amplitudes{run_number} = run_amplitudes;
end
if plot_results
    figure
end
hold on
for run_number =run_numbers
    plot(sensor_locations, amplitudes{run_number}, '-x','MarkerSize', 20,'DisplayName',sprintf("%s run number %d", ex_day, run_number))
end
xlabel('sensor distance from piston [m]')
ylabel("measured amplitude [m]")
title('measured amplitudes different sensors')
fontsize(20, "points");
legend()

%% run number 1
lab_run_number = 1;
basilisk_folders = [];
titles = [];
runtimes = [];

%% run 1
run_numbers = [1 2 3];
basilisk_folders =[];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL11_layers10/"]; titles = [titles "multilayer LEVEL:11, nl:10"]; runtimes=[runtimes 843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL12_layers10/"]; titles = [titles "multilayer LEVEL:12, nl:10"]; runtimes=[runtimes 843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL11/"]; titles = [titles "boundary LEVEL 11"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL12/"]; titles = [titles "boundary LEVEL 12"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/fil2/LEVEL11/"]; titles = [titles "fil2 boundary LEVEL 11"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/fil2/LEVEL12/"]; titles = [titles "fil2 boundary LEVEL 12"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL11_0/"]; titles = [titles "moving piston LEVEL 11"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL12_0/"]; titles = [titles "moving piston LEVEL 12"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL12_1/"]; titles = [titles "moving piston LEVEL 12+1"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL12_2/"]; titles = [titles "moving piston LEVEL 12+2"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL13_0/"]; titles = [titles "moving piston LEVEL 13"];runtimes=[runtimes 6843];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL13_1/"]; titles = [titles "moving piston LEVEL 13+1"];runtimes=[runtimes 6843];



%% run 4
run_numbers = [6];
basilisk_folders =[];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run4/LEVEL11_layers10/"]; titles = [titles "multilayer LEVEL:11, nl:10"]; runtimes=[runtimes 843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run4/LEVEL12_layers10/"]; titles = [titles "multilayer LEVEL:12, nl:20"]; runtimes=[runtimes 843];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run4/LEVEL13_layers20/"]; titles = [titles "multilayer LEVEL:13, nl:20"]; runtimes=[runtimes 843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run4/LEVEL11/"]; titles = [titles "boundary LEVEL 11"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run4/LEVEL12/"]; titles = [titles "boundary LEVEL 12"];runtimes=[runtimes 6843];
%basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run4/LEVEL13/"]; titles = [titles "boundary LEVEL 13"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run4/LEVEL11_0/"]; titles = [titles "moving piston LEVEL 11"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run4/LEVEL12_0/"]; titles = [titles "moving piston LEVEL 12"];runtimes=[runtimes 6843];
%% run 5
run_numbers = [7 8 9];
basilisk_folders =[];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run5/LEVEL11_layers10/"]; titles = [titles "mulitlayer LEVEL:11, nl:10"]; runtimes=[runtimes 843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run5/LEVEL12_layers10/"]; titles = [titles "multilayer LEVEL:12, nl:10"]; runtimes=[runtimes 843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run5/LEVEL11/"]; titles = [titles "boundary LEVEL 11"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run5/LEVEL12/"]; titles = [titles "boundary LEVEL 12"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run5/LEVEL13/"]; titles = [titles "boundary LEVEL 13"];runtimes=[runtimes 6843];

basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run5/LEVEL11_0/"]; titles = [titles "moving piston LEVEL 11"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run5/LEVEL12_0/"]; titles = [titles "moving piston LEVEL 12"];runtimes=[runtimes 6843];

%%
total_probes = 1404;
basilisk_start_indices = [start_indices(lab_run_number,:) start_indices(lab_run_number,1):floor((start_indices(lab_run_number,end)-start_indices(lab_run_number,1))/(total_probes-5)):start_indices(lab_run_number,end)];
basilisk_start_indices = [start_indices(lab_run_number,:) 2000:floor((4000-2000)/(total_probes-5)):4000];

first_probe = 21; %probes for plotting and finding amplitudes
last_probe = 1404;
basilisk_probes = first_probe:last_probe;
% start_indices(lab_run_number, 1) = 1;
% basilisk_start_indices = [zeros(1,basilisk_probes(1)-1) start_indices(lab_run_number,1):floor((start_indices(lab_run_number,end)-start_indices(lab_run_number,1))/(size(basilisk_probes,2)-1)):start_indices(lab_run_number,end)];
% basilisk_start_indices = basilisk_start_indices*0+1;
basilisk_end_indices = min(basilisk_start_indices+3600, 5000);
%%

close all
%plot the basilisk results
Fs = 100;
f = 1.425;
min_separation = 1/f*Fs*0.8; %minimum separation when finding the crests and troughs
plot_results = false;
basilisk_wave_amplitudes = {};
for i=1:size(basilisk_folders, 2)
% for i=1
    
    basilisk_folder = basilisk_folders(i);
    %plot the basilisk results
    surface_probes = readtable(append(basilisk_folder, "surface_probes.csv"), 'ReadVariableNames',true, 'VariableNamingRule','preserve');
    probe_names = surface_probes.Properties.VariableNames;
    probe_names{20} = probe_names{20}(1:8);
    probe_names{120} = probe_names{120}(1:8);
    basilisk_probe_locations = str2double(probe_names');
    basilisk_probe_locations = basilisk_probe_locations(2:end);
    surface_probes = table2array(surface_probes);
    if false
        % using basilisk gauges instead
        x_start = 1.5;
        x_end = 14;
        gauge_mat = [];
        for x=x_start:0.1:x_end
            gauge = load(sprintf("%sx%.3f", folder, x));
            gauge_mat = [gauge_mat gauge(:,2)];
            plot(gauge(:,1)*100, gauge(:,2), '--', 'DisplayName', sprintf("%.2f", x))
        end
        surface_probes = gauge_mat;
    end
    if plot_results
        figure
        hold on 
        for sensor = basilisk_probes
            % plot(surface_probes(basilisk_start_indices(sensor):basilisk_end_indices(sensor),sensor+1), DisplayName=sprintf("probe at x=%.1fm", basilisk_probe_locations(sensor)))
            plot(surface_probes(:,sensor+1), DisplayName=sprintf("probe at x=%.1fm", basilisk_probe_locations(sensor)))
        end
        title(titles(i));
        legend()
    end
    basilisk_run_amplitudes = [];
    if plot_results
        figure;
        hold on
    end
    for sensor=basilisk_probes
        basilisk_probe_number = sensor+1;
        basilisk_surface = surface_probes(basilisk_start_indices(sensor):basilisk_end_indices(sensor),basilisk_probe_number);
        t = 0:1/Fs:(size(surface_probes(:,basilisk_probe_number), 1)-1)/Fs;
        t = t(basilisk_start_indices(sensor):basilisk_end_indices(sensor));
        LMax = islocalmax(basilisk_surface, 1, 'MinSeparation',min_separation);
        LMin = islocalmin(basilisk_surface, 1, 'MinSeparation',min_separation);
        n_maxima = size(basilisk_surface(LMax),1);
        n_minima = size(basilisk_surface(LMin),1);
        if n_maxima>n_minima
            LMax(find(LMax,1, 'last')) = 0;
        elseif n_maxima<n_minima
            LMin(find(LMin,1, 'last')) = 0;
        end
        n_maxima = size(basilisk_surface(LMax),1);
        n_minima = size(basilisk_surface(LMin),1);
        if n_maxima>n_minima
            LMax(find(LMax,1, 'first')) = 0;
        elseif n_maxima<n_minima
            LMin(find(LMin,1, 'first')) = 0;
        end
        
        basilisk_run_amplitude = (basilisk_surface(LMax)-basilisk_surface(LMin))/2;
        asdf = mean(basilisk_run_amplitude);
        basilisk_run_amplitudes = [basilisk_run_amplitudes mean(basilisk_run_amplitude)];
        if plot_results
            if true
                plot(t(LMax), basilisk_run_amplitude,'x-', "DisplayName",sprintf("probe %d", sensor))
            else
                plot(t, basilisk_surface)
                plot(t(LMax), basilisk_surface(LMax), 'x')
                plot(t(LMin), basilisk_surface(LMin), 'x')
            end
        end
    end
    if plot_results
        title(titles(i))
        legend()
    end
    basilisk_wave_amplitudes{i} = basilisk_run_amplitudes;
end
%plot the evolution of the amplitude
% close all
figure
hold on
for i=1:size(basilisk_folders,2)
    plot(basilisk_probe_locations(basilisk_probes)', basilisk_wave_amplitudes{i}, 'DisplayName',titles(i))
end
for run_number=run_numbers
    plot(sensor_locations, amplitudes{run_number}, '-x','MarkerSize', 15, 'DisplayName',sprintf("Lab experiment run number %d", run_number))
end
fontsize(20, "points")
title("Mean amplitude at distance from piston")
xlabel("Distance piston to sensor [m]")
ylabel("mean amplitude at sensor position [m]")
legend()


%%
folder = "~/Documents/master/basilisk/2d_piston/multilayer-piston/";
% figure 
hold on
x_start = 1.5;
x_end = 14;
gauge_mat = [];
for x=x_start:0.1:x_end
    gauge = load(sprintf("%sx%.3f", folder, x));
    gauge_mat = [gauge_mat gauge(:,2)];
    plot(gauge(:,1)*100, gauge(:,2), '--', 'DisplayName', sprintf("%.2f", x))
end
legend


%%


h = 0.6;
a = 0.01;
T = 1/1.425;
Result = StokesDispSolver('h', h, 'H', 2*a, 'T', T, 'mode', 1);
theta = 0:0.01:40*2*pi;
[eta, eta1, eta2, eta3,eta4,eta5] = StokesEta(Result.k, h, Result.a, theta);
LMax= islocalmax(eta);
LMin = islocalmin(eta);
n_maxima = size(eta(LMax),2);
n_minima = size(eta(LMin),2);
if n_maxima>n_minima
    LMax(find(LMax,1, 'last')) = 0;
elseif n_maxima<n_minima
    LMin(find(LMin,1, 'last')) = 0;
end

t_sensor = {};


wave_amplitudes = {};

wave_amplitudes{1} = (eta(LMax)-eta(LMin))/2;
plot(theta(LMax)/2/pi, wave_amplitudes{1})
