clear all
close all
lab_files
%plot with basilisk results
%% run number 7
lab_run_number = 7;
basilisk_folders = [];
titles = [];
runtimes = [];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL11_0/"]; titles = [titles "LEVEL 11"]; runtimes=[runtimes 843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run5/LEVEL11/"]; titles = [titles "boundary LEVEL 11"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run5/LEVEL12/"]; titles = [titles "boundary LEVEL 12"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run5/LEVEL13/"]; titles = [titles "boundary LEVEL 13"];runtimes=[runtimes 6843];

%% run number 1
lab_run_number = 1;
basilisk_folders = [];
titles = [];
runtimes = [];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL11_0/"]; titles = [titles "LEVEL 11"]; runtimes=[runtimes 843];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL11_1/"]; titles = [titles "LEVEL 11+1"];runtimes=[runtimes 1792];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL11_2/"]; titles = [titles "LEVEL 11+2"];runtimes=[runtimes 4308];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL12_0/"]; titles = [titles "LEVEL 12+0"];runtimes=[runtimes 4089];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL12_1/"]; titles = [titles "LEVEL 12+1"];runtimes=[runtimes 4089];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL13_0/"]; titles = [titles "LEVEL 13+0"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/no_DT/LEVEL11/"]; titles = [titles "boundary LEVEL 11"];runtimes=[runtimes 350];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/no_DT/LEVEL12/"]; titles = [titles "boundary LEVEL 12"];runtimes=[runtimes 6843];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/no_DT/LEVEL13/"]; titles = [titles "boundary LEVEL 13"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/conserving/LEVEL11/"]; titles = [titles "conserving LEVEL 11"];runtimes=[runtimes 350];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/conserving/LEVEL12/"]; titles = [titles "conserving LEVEL 12"];runtimes=[runtimes 6843];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/conserving/LEVEL13/"]; titles = [titles "conserving LEVEL 13"];runtimes=[runtimes 6843];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/no_conserving_with_DT/LEVEL11/"]; titles = [titles "DT l/2**LEVEL, LEVEL 11"];runtimes=[runtimes 350];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL11/"]; titles = [titles "conserving.h DT l/2**LEVEL/2, LEVEL 11"];runtimes=[runtimes 350];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL12/"]; titles = [titles "conserving.h DT l/2**LEVEL/2, LEVEL 12"];runtimes=[runtimes 350];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL13/"]; titles = [titles "conserving.h DT l/2**LEVEL/2, LEVEL 13"];runtimes=[runtimes 350];
%% plot one sensor from experiment with the basilisk results
sensor= 1;
basilisk_probe_number = sensor+1;
figure;
hold on
lab_run_number=3
for lab_run_number=lab_run_number
    file = load(append(lab_folder, lab_surface_probes(lab_run_number)));
    t = 0.00:0.008:(size(file, 1)-1)/125;



    %plot the lab results
    col = sensor+2;
    plot(t, -(file(:,col)-file(1,col)), '.', 'DisplayName',sprintf('run %d surface probe', lab_run_number))
end

for i=1:size(basilisk_folders, 2)
    basilisk_folder = basilisk_folders(i);
    %plot the basilisk results
    
    % surface_probes = load(append(basilisk_folder, "surface_probes.csv"));
    surface_probes = readtable(append(basilisk_folder, "surface_probes.csv"), 'ReadVariableNames',true, 'VariableNamingRule','preserve');
    probe_locations = str2double(surface_probes.Properties.VariableNames);
    surface_probes = table2array(surface_probes);
    plot(surface_probes(:,1), surface_probes(:,basilisk_probe_number),'DisplayName',titles(i));
    xlabel('t[s]')
    ylabel('surface eleveation[m]')
    % xlim([0 10])
    %mark the period
    %plot(0:T:60,0, 'x')
end
legend();
title(sprintf("surface elevation, sensor %.2fm from piston at rest", probe_locations(basilisk_probe_number)))

%% plot sensor 1 for different basilisk runs in different plots
sensor = 1;
file = load(append(lab_folder, lab_surface_probes(lab_run_number)));
t = 0.008:0.008:size(file, 1)/125;

figure;
tiledlayout(size(basilisk_folders, 2),1)
xmin = 1;
xmax = 9;
for i=1:size(basilisk_folders, 2)
    basilisk_folder = basilisk_folders(i);
    titl = titles(i);
    nexttile;
    hold on

    %plot the lab results
    col = sensor+2;
    max(file(:,col)) - min(file(:,col))/2;
    plot(t, -(file(:,col)-file(1,col)), 'DisplayName','surface probe')
    
    %plot the basilisk results
    surface_probes = readtable(append(basilisk_folder, "surface_probes.csv"));
    surface_probes = table2array(surface_probes);
    plot(surface_probes(:,1), surface_probes(:,sensor+1),'DisplayName','basilisk')
    title(titl)
    legend()
    xlabel('t[s]')
    ylabel('surface eleveation[m]')
    xlim([xmin xmax])
    %mark the period
    %plot(0:T:60,0, 'x')
end
%% plot all the sensors for one run
%basilisk_folder = "~/Documents/master/basilisk_results/a_308/LEVEL13_no_tension/piston-moving/"; lab_run_number=1;
%basilisk_folder = "~/Documents/master/basilisk/piston-moving/"; lab_run_number=5;
run_number = 2;
basilisk_folder = basilisk_folders(run_number);


omega = 8.95; 
T = (2*pi)/omega;
%load the lab results
file = load(append(lab_folder, lab_surface_probes(lab_run_number)));
t = 0.008:0.008:size(file, 1)/125;
probe_positions = [8 10 10.7 11.5];
%plot the lab and basilisk results
surface_probes = load(append(basilisk_folder, "surface_probes.csv"));
figure;
tiledlayout(4,1);

for sensor=1:4
    nexttile;
    hold on
    
    %plot the lab results
    col = sensor+2;
    max(file(:,col)) - min(file(:,col))/2;
    plot(t, -(file(:,col)-file(1,col)), 'DisplayName','surface probe')

    %plot the basilisk results
    plot(surface_probes(:,1), surface_probes(:,sensor+1),'DisplayName','basilisk')
    title(sprintf("sensor %d, x=%.1f", sensor, probe_positions(sensor)))
    legend()
    xlabel('t[s]')
    ylabel('surface eleveation[m]')
    xlim([0 30])
    %mark the period
    %plot(0:T:60,0, 'x')
end


%%
figure;
hold on;
for lab_run_number = 1:5
    load(sprintf("Lab_results/24_04_12/%d/padle_ut.dat", lab_run_number))
    plot((padle_ut-padle_ut(1))*4.2, 'DisplayName',sprintf("padle ut run number:%d",lab_run_number))
    load(sprintf("Lab_results/24_04_12/%d/fil3.dat", lab_run_number))
    plot(fil3-fil3(1),'--', 'DisplayName',sprintf("fil3 run number:%d",lab_run_number))
end
legend()
title("padle_ut")
%%
figure;
hold on;
for lab_run_number = 1
    load(sprintf("Lab_results/24_09_18/%d/fil1.dat", lab_run_number))
    plot(-(fil1-fil1(1))*4.4, 'DisplayName',sprintf("fil1 run number:%d",lab_run_number))
    load(sprintf("Lab_results/24_04_12/%d/fil1.dat", lab_run_number))
    plot(-(fil1-fil1(1))*4.4, 'DisplayName',sprintf("fil1 run number APRIL:%d",lab_run_number))
    % load(sprintf("Lab_results/24_04_12/%d/fil3.dat", lab_run_number))
    % plot(fil3-fil3(1),'--', 'DisplayName',sprintf("fil3 run number:%d",lab_run_number))
    load(sprintf("Lab_results/24_04_12/%d/piston_position.dat", lab_run_number))
    plot((piston_position-piston_position(1))*100,'-', 'DisplayName',sprintf("piston_position smoothed run number:%d",lab_run_number))
    load(sprintf("Lab_results/24_04_12/%d/piston_speed.dat", lab_run_number))
    plot((piston_speed)*100,'-', 'DisplayName',sprintf("piston_speed run number:%d",lab_run_number))
end
xlim([0, 600])
legend()
title("fil1")
%% plot fil1 fil3 and piston_position

figure
for lab_run_number = 3
    figure;
    hold on;
    load(sprintf("Lab_results/24_04_12/%d/fil1.dat", lab_run_number))
    load(sprintf("Lab_results/24_04_12/%d/fil2.dat", lab_run_number))
    load(sprintf("Lab_results/24_04_12/%d/fil3.dat", lab_run_number))
    load(sprintf("Lab_results/24_04_12/%d/padle_ut.dat", lab_run_number))
    load(sprintf("Lab_results/24_04_12/%d/piston_position.dat", lab_run_number))

    %load(sprintf("Lab_results/24_09_18/%d/fil1.dat", lab_run_number))
    %load(sprintf("Lab_results/24_09_18/%d/fil3.dat", lab_run_number))
    %load(sprintf("Lab_results/24_09_18/%d/piston_position.dat", lab_run_number))

    fil1 = -(fil1 - fil1(1));
    plot(fil1*0.044, 'DisplayName',sprintf("fil1 run number:%d",lab_run_number))
    
    fil2 = fil2 - fil2(1);
    % plot(fil2*0.01,'-', 'DisplayName',sprintf("fil2 run number:%d",lab_run_number))

    fil3 = (fil3-fil3(1));
    plot(fil3*0.01,'-', 'DisplayName',sprintf("fil3 run number:%d",lab_run_number))

    padle_ut = padle_ut-padle_ut(1);
    plot(padle_ut*0.044,'-', 'DisplayName',sprintf("padle ut run number:%d",lab_run_number))

    plot(piston_position,'-', 'DisplayName',sprintf("piston_position run number:%d",lab_run_number))
    legend()
    title(sprintf("run %d", lab_run_number))

end

%%
figure;
hold on;
for lab_run_number = 1:5
    load(sprintf("Lab_results/24_04_12/%d/fil2.dat", lab_run_number))
    plot(fil2-fil2(1), 'DisplayName',sprintf("run number:%d",lab_run_number))
end
legend()
title("fil2")

figure;
hold on;
for lab_run_number = 1:5
    load(sprintf("Lab_results/24_04_12/%d/fil3.dat", lab_run_number))
    plot(fil3-fil3(1), 'DisplayName',sprintf("run number:%d",lab_run_number))
end
legend()
title("fil3")
