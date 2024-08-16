clear all
close all
%% plot the piston movements
lab_folder_number = 1;
lab_folder = sprintf("Lab_results/24_04_12/%d/" , lab_folder_number);
figure;
hold on;
load(append(lab_folder, "fil1.dat"));
fil = fil1;
t = 0.008:0.01:size(fil1, 1)/100;

plot(t, fil-fil(1,1))
load(append(lab_folder, "fil2.dat"));
fil = fil2;
plot(t, fil-fil(1,1))

load(append(lab_folder, "fil3.dat"));
fil = fil3;
plot(t, fil-fil(1,1))

load(append(lab_folder, "padle_ut.dat"));
plot(t, (padle_ut-padle_ut(1)))

legend('fil1', 'fil2', 'fil3', 'padle ut')

hs = 0:1/1.425:60;
xs = zeros(size(hs));
plot(hs,xs, 'X');
title("plot of the paddle files");

%% plot with basilisk results
close all
%% results for the wave from run number 5
%half amplitude piston movement
lab_run_number = 5;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL12/"]; titles = [titles "LEVEL 12"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL13_2pad/"]; titles = [titles "LEVEL 13 2 padding"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL13_refined_bottom/"]; titles = [titles "LEVEL 13 refined bottom"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL13_unrefined_bottom/"]; titles = [titles "LEVEL 13 unrefined bottom"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14/"]; titles = [titles "LEVEL 14"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_1pad/"]; titles = [titles "LEVEL 14 1 pad"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL15/"]; titles = [titles "LEVEL 15"];

%% wave 5 testing refinement parameters
%half amplitude piston movement
lab_run_number = 5;
basilisk_folders = [];
titles = [];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL13_refineparam1/"]; titles = [titles "LEVEL 13 refineparam 1"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL13_refineparam01/"]; titles = [titles "LEVEL 13 refineparam 0.1"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL13_refineparam001/"]; titles = [titles "LEVEL 13 refineparam 0.01"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL13_refineparam0001/"]; titles = [titles "LEVEL 13 refineparam 0.001"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14/"]; titles = [titles "LEVEL 14 refineparam 1"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam1/"]; titles = [titles "LEVEL 14 refineparam 1"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam01/"]; titles = [titles "LEVEL 14 refineparam 0.1"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam001/"]; titles = [titles "LEVEL 14 refineparam 0.01"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam0001/"]; titles = [titles "LEVEL 14 refineparam 0.001"];
%% run 4 wave results
%double amplitude piston movement
lab_run_number = 4;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL12/"]; titles = [titles "LEVEL 12"];
basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL13/"]; titles = [titles "LEVEL 13"];
basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL13_2pad/"]; titles = [titles "LEVEL 13 2 padding"];
basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL14_1pad/"]; titles = [titles "LEVEL 14 1 padding"];
basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL14_tension/"]; titles = [titles "LEVEL 14 with surface tension"];
basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL_14/"]; titles = [titles "LEVEL 14"];

%% run 1 wave results
%standard piston amplitude
lab_run_number = 1;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/results/run1/LEVEL12/"]; titles = [titles "LEVEL 12"];
basilisk_folders = [basilisk_folders "~/Documents/results/run1/LEVEL13/"]; titles = [titles "LEVEL 13"];
basilisk_folders = [basilisk_folders "~/Documents/results/run1/LEVEL13_2pad/"]; titles = [titles "LEVEL 13 2 padding"];

%% plot all the sensors for one run
%basilisk_folder = "~/Documents/master/basilisk_results/a_308/LEVEL13_no_tension/piston-moving/"; lab_run_number=1;
%basilisk_folder = "~/Documents/master/basilisk/piston-moving/"; lab_run_number=5;
run_number = 1;
basilisk_folder = basilisk_folders(run_number);


omega = 8.95; 
T = (2*pi)/omega;
%load the lab results
lab_folder = "Lab_results/24_04_12/";
filenames = ["f1425_a0_308.csv";
            "f1425_a0_308_r2.csv";
            "f1425_a0_308_r3.csv";
            "f1425_a0_616_r1.csv";
            "f1425_a0_154_r1.csv"];
file = load(append(lab_folder, filenames(lab_run_number)));
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
    plot(t, file(:,col)-file(1,col), 'DisplayName','surface probe')

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

%% plot sensor 1 for different basilisk runs in different plots
sensor = 1;
figure;
tiledlayout(size(basilisk_folders, 2),1)

for i=1:size(basilisk_folders, 2)
    basilisk_folder = basilisk_folders(i);
    titl = titles(i);
    nexttile;
    hold on

    %plot the lab results
    col = sensor+2;
    max(file(:,col)) - min(file(:,col))/2;
    plot(t, file(:,col)-file(1,col), 'DisplayName','surface probe')
    
    %plot the basilisk results
    surface_probes = load(append(basilisk_folder, "surface_probes.csv"));
    plot(surface_probes(:,1), surface_probes(:,sensor+1),'DisplayName','basilisk')
    title(titl)
    legend()
    xlabel('t[s]')
    ylabel('surface eleveation[m]')
    xlim([0 30])
    %mark the period
    %plot(0:T:60,0, 'x')
end
%% plot sensor 1 for the results in basilisk_folders in the same plot
sensor=1;
figure;
hold on

%plot the lab results
col = sensor+2;
max(file(:,col)) - min(file(:,col))/2;
plot(t, file(:,col)-file(1,col), '--', 'DisplayName','surface probe')


for i=1:size(basilisk_folders, 2)
    basilisk_folder = basilisk_folders(i);
    %plot the basilisk results
    surface_probes = load(append(basilisk_folder, "surface_probes.csv"));
    plot(surface_probes(:,1), surface_probes(:,sensor+1),'DisplayName',titles(i));
    xlabel('t[s]')
    ylabel('surface eleveation[m]')
    xlim([0 30])
    %mark the period
    %plot(0:T:60,0, 'x')
end
legend();

%%
figure;
hold on;
for lab_run_number = 1:5
    load(sprintf("Lab_results/24_04_12/%d/padle_ut.dat", lab_run_number))
    plot(padle_ut-padle_ut(1), 'DisplayName',sprintf("run number:%d",lab_run_number))
end
legend()
title("padle_ut")

figure;
hold on;
for lab_run_number = 1:5
    load(sprintf("Lab_results/24_04_12/%d/fil1.dat", lab_run_number))
    plot(fil1-fil1(1), 'DisplayName',sprintf("run number:%d",lab_run_number))
end
legend()
title("fil1")

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
