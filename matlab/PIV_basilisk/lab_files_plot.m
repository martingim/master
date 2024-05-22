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

%basilisk_folder ="~/2piston/"; lab_run_number = 4;
basilisk_folder = "~/piston-moving/"; lab_run_number=2;

%load the lab results
lab_folder = "Lab_results/24_04_12/";
filenames = ["f1425_a0_308.csv";
            "f1425_a0_308_r2.csv";
            "f1425_a0_308_r3.csv";
            "f1425_a0_616_r1.csv";
            "f1425_a0_154_r1.csv"];
file = load(append(lab_folder, filenames(lab_run_number)));
t = 0.008:0.008:size(file, 1)/125;
%plot the lab results
figure;
hold on

for col=3:3
    max(file(:,col)) - min(file(:,col))/2;
    plot(t, file(:,col)-file(1,col))
end

legend('sensor1', 'sensor2', 'sensor3', 'sensor4', 'sensor5', 'sensor6')
% 
% load basilisk_results/X_0
% plot(X_0(:,1), X_0(:,2), 'DisplayName','local basilisk')
% 
% load basilisk_results/X_0_1
% plot(X_0_1(:,1), X_0_1(:,2), 'DisplayName','abacus basilisk')
% 

%load and plot the basilsk results
surface_probes = load(append(basilisk_folder, "surface_probes.csv"));
plot(surface_probes(:,1), surface_probes(:,2),'DisplayName','NS')
title(sprintf("lab run number:%d", lab_run_number));

% X marks the period 
% hs = 0:1/1.425:60;
% xs = zeros(size(hs));
% plot(hs,xs, 'X');

%%
load Lab_results/24_04_12/4/padle_ut.dat
load Lab_results/24_04_12/4/fil3.dat
%plot((padle_ut-padle_ut(1)))

plot(fil3-fil3(1))