clear all
close all
lab_files
%plot with basilisk results
%% results for the wave from run number 5
%half amplitude piston movement
lab_run_number = 7;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL12/"]; titles = [titles "LEVEL 12"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL13/"]; titles = [titles "LEVEL 13"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14/"]; titles = [titles "LEVEL 14"];
%% results for the wave from run number 5 test LEVEL
%half amplitude piston movement
lab_run_number = 7;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL12/"]; titles = [titles "LEVEL 12"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/ LEVEL13_refineparam0001/"]; titles = [titles "LEVEL 13"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam0001/"]; titles = [titles "LEVEL 14"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL15/"]; titles = [titles "LEVEL 15"];

%% wave 5 testing refinement parameters
%half amplitude piston movement
lab_run_number = 7;
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
% basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL12/"]; titles = [titles "LEVEL 12"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL13/"]; titles = [titles "LEVEL 13"];
basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL14/"]; titles = [titles "LEVEL 14"];

%% run 1 wave results
%standard piston amplitude
lab_run_number = 1;
basilisk_folders = [];
titles = [];
% basilisk_folders = [basilisk_folders "~/Documents/results/run1/LEVEL12/"]; titles = [titles "LEVEL 12"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run1/LEVEL13/"]; titles = [titles "LEVEL 13"];
basilisk_folders = [basilisk_folders "~/Documents/results/run1/LEVEL13/"]; titles = [titles "LEVEL 13"];
basilisk_folders = [basilisk_folders "~/Documents/results/run1/LEVEL13_15/"]; titles = [titles "LEVEL 13 piston 15"];

%% run 1 compare boundary to piston set 
%standard piston amplitude
lab_run_number = 1;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/LEVEL10_0/"]; titles = [titles "LEVEL 10 no adaptive refinement"];
basilisk_folders = [basilisk_folders "~/Documents/LEVEL10_adaptive/"]; titles = [titles "LEVEL 10 adaptive refinement"];
basilisk_folders = [basilisk_folders "~/Documents/boundary/"]; titles = [titles "LEVEL 10 set boundary velocity"];

%% multilayer compare ratio of height to width of the cells in the top layer
%standard piston amplitude
lab_run_number = 1;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL10_layers10/"]; titles = [titles "multilayer LEVEL 10 nl=10"];
%basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL11_layers5/"]; titles = [titles "multilayer LEVEL 11 height/width=3.44"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL11_layers10/"]; titles = [titles "multilayer LEVEL 11 nl=10"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL11_layers20/"]; titles = [titles "multilayer LEVEL 11 nl=20"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL12_layers20/"]; titles = [titles "multilayer LEVEL 12 nl=20"];


%% compare multilayer to NS for same LEVEL
%standard piston amplitude
lab_run_number = 1;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL9_layers10/"]; titles = [titles "multilayer LEVEL 9 nl=10"];
%basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL11_layers5/"]; titles = [titles "multilayer LEVEL 11 height/width=3.44"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL11_layers10/"]; titles = [titles "multilayer LEVEL 11 nl=10"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL11_layers20/"]; titles = [titles "multilayer LEVEL 11 nl=20"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL10/"]; titles = [titles "NS LEVEL 10"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL9/"]; titles = [titles "NS LEVEL 9"];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL11/"]; titles = [titles "NS LEVEL 11"];

%% run number 2, 2d boundary piston
lab_run_number = 2;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run2/LEVEL10/"]; titles = [titles "NS LEVEL 10"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run2/LEVEL9/"]; titles = [titles "NS LEVEL 9"];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL11/"]; titles = [titles "NS LEVEL 11"];


%% run number 3, 2d boundary piston
lab_run_number = 3;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL10/"]; titles = [titles "run 1 LEVEL 10"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run2/LEVEL10/"]; titles = [titles "run 2 LEVEL 10"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run3/LEVEL10/"]; titles = [titles "run 3 LEVEL 10"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run3/LEVEL11/"]; titles = [titles "run 3 LEVEL 11"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run3/LEVEL12/"]; titles = [titles "run 3 LEVEL 12"];

%% run number 3, 2d boundary piston
lab_run_number = 4;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run4/LEVEL10/"]; titles = [titles "NS LEVEL 10"];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run3/LEVEL11/"]; titles = [titles "NS LEVEL 11"];

%% run number 7, 2d boundary piston
lab_run_number = 7;
basilisk_folders = [];
titles = [];

basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run7/LEVEL10/"]; titles = [titles "NS LEVEL 10"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run7/LEVEL11/"]; titles = [titles "NS LEVEL 11"];

%% run number 1
lab_run_number = 1;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL10/"]; titles = [titles "NS LEVEL 10"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL11_0/"]; titles = [titles "moving piston LEVEL 11"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/piston-moving/results/run1/LEVEL12_0/"]; titles = [titles "moving piston LEVEL 12"];

%basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run3/LEVEL10/"]; titles = [titles "NS LEVEL 10"];

%% plot sensor 1 for different basilisk runs in different plots
sensor = 2;
file = load(append(lab_folder, lab_surface_probes(lab_run_number)));
t = 0.008:0.008:size(file, 1)/125;

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
    plot(t, -(file(:,col)-file(1,col)), 'DisplayName','surface probe')
    
    %plot the basilisk results
    surface_probes = readtable(append(basilisk_folder, "surface_probes.csv"));
    surface_probes = table2array(surface_probes);
    plot(surface_probes(:,1), surface_probes(:,sensor+1),'DisplayName','basilisk')
    title(titl)
    legend()
    xlabel('t[s]')
    ylabel('surface eleveation[m]')
    xlim([22 28])
    %mark the period
    %plot(0:T:60,0, 'x')
end
%% plot one sensor from experiment with the basilisk results
sensor=1; 
basilisk_probe_number = sensor+1;
figure;
hold on
for lab_run_number=lab_run_number
    file = load(append(lab_folder, lab_surface_probes(lab_run_number)));
    t = 0.00:0.008:(size(file, 1)-1)/125;



    %plot the lab results
    col = sensor+2;
    plot(t, -(file(:,col)-file(1,col)), 'x', 'DisplayName',sprintf('run %d surface probe', lab_run_number))
end

for i=1:size(basilisk_folders, 2)
    basilisk_folder = basilisk_folders(i);
    %plot the basilisk results
    
    % surface_probes = load(append(basilisk_folder, "surface_probes.csv"));
    surface_probes = readtable(append(basilisk_folder, "surface_probes.csv"));
    surface_probes = table2array(surface_probes);
    plot(surface_probes(:,1), surface_probes(:,basilisk_probe_number),'DisplayName',titles(i));
    xlabel('t[s]')
    ylabel('surface eleveation[m]')
    xlim([0 10])
    %mark the period
    %plot(0:T:60,0, 'x')
end
legend();
title(sprintf("surface elevation for sensor %d", sensor))


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
    plot(fil2*0.01,'-', 'DisplayName',sprintf("fil2 run number:%d",lab_run_number))

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
