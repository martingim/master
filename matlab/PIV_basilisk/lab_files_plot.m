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
basilisk_folders = [];
titles = [];
%half amplitude piston movement
%basilisk_folder = "~/Documents/master/basilisk_results/a_154/LEVEL15_embed_reduced/"; lab_run_number=5;
%basilisk_folder = "~/Documents/master/basilisk_results/a_154/LEVEL13_embed_reduced/"; lab_run_number=5;
%basilisk_folder = "~/Documents/master/basilisk_results/a_154/LEVEL14/"; lab_run_number=5;

%308 piston amplitude
lab_run_number = 1;
%basilisk_folder = "~/Documents/master/basilisk_results/a_308/LEVEL/"; lab_run_number=2;
%normal amplitude piston movement
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk_results/a_308/LEVEL10/"]; titles = [titles "LEVEL 10, 2.4 cm cells"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk_results/a_308/LEVEL11/"]; titles = [titles "LEVEL 11, 1.2 cm cells"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk_results/a_308/LEVEL12/"]; titles = [titles "LEVEL 12, 0.6 cm cells"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk_results/a_308/LEVEL13/"]; titles = [titles "LEVEL 13, 0.3 cm cells"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk_results/a_308/LEVEL13_2/"]; titles = [titles "LEVEL 13_2"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk_results/a_308/LEVEL13_wavetank/"]; titles = [titles "LEVEL 13 not reduced"];

basilisk_folder = "~/Documents/master/basilisk_results/a_308/LEVEL13_no_tension/piston-moving/"; lab_run_number=1;


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

%% plot sensor 1 for different basilisk runs
sensor = 1;
figure;
tiledlayout(6,1)

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
