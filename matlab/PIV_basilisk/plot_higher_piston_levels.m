lab_files
%%
lab_run_number = 1;
run_number = 1;
basilisk_folders = [];
titles = [];
% basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL10_1/"]; titles = [titles "LEVEL 10 Piston + 1, 0:00:55"];
% basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL10_2/"]; titles = [titles "LEVEL 10 Piston + 2, 0:02:09"];
% basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL11/"]; titles = [titles "LEVEL 11 Piston + 0, 0:01:48"];
% basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL11_1/"]; titles = [titles "LEVEL 11 Piston + 1"];
% basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL11_2/"]; titles = [titles "LEVEL 11 Piston + 2, 0:10:35"];
% basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL11_3/"]; titles = [titles "LEVEL 11 Piston + 3"];
% basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL12/"]; titles = [titles "LEVEL 12 Piston + 0, 0:10:27"];
% basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL12_1/"]; titles = [titles "LEVEL 12 Piston + 1, 0:29:26"];
% basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL12_2/"]; titles = [titles "LEVEL 12 Piston + 2, 1:19:30"];
% basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL13/"]; titles = [titles "LEVEL 13 Piston + 0, 0:57:23"];
% basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL13_1/"]; titles = [titles "LEVEL 13 Piston + 1, 3:46:40"];
basilisk_folders = [basilisk_folders "~/Documents/results/test_two_levels/LEVEL14/"]; titles = [titles "NS LEVEL 14 Piston set water velocity"];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/boundary-piston/results/run1/LEVEL10/"]; titles = [titles "NS LEVEL 10 set boundary"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/boundary-piston/results/run1/LEVEL11/"]; titles = [titles "NS LEVEL 11 set boundary"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/boundary-piston/results/run1/LEVEL12/"]; titles = [titles "LEVEL 12 set boundary"]; 
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/multilayer-piston/results/run1/LEVEL10_layers10/"]; titles = [titles "multilayer LEVEL 10 set boundary"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/multilayer-piston/results/run1/LEVEL11_layers20/"]; titles = [titles "multilayer LEVEL 11 set boundary"];
% basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/multilayer-piston/results/run1/LEVEL12_layers40/"]; titles = [titles "multilayer LEVEL 12 set boundary"];


file = load(append(lab_folder, lab_surface_probes(lab_run_number)));
t = 0.008:0.008:size(file, 1)/125;
% t = 0.008:0.01:size(fil1, 1)/100;
%% plot sensor 1 at 1.5m for the results in basilisk_folders in the same plot
file = load(append(lab_folder, lab_surface_probes(lab_run_number)));
sensor=1;
basilisk_probe_number = 21;

figure;
hold on
title("Surface probe at 1.5m")
%plot the lab results
col = sensor+2;
max(file(:,col)) - min(file(:,col))/2;
plot(t, -(file(:,col)-file(1,col)), '--', 'DisplayName','surface probe wavetank')
i = 1;

for i=1:size(basilisk_folders, 2)
    basilisk_folder = basilisk_folders(i);
    %plot the basilisk results
    
    % surface_probes = load(append(basilisk_folder, "surface_probes.csv"));
    surface_probes = readtable(append(basilisk_folder, "surface_probes.csv"), ReadVariableNames=true);
    surface_probes = table2array(surface_probes);
    plot(surface_probes(:,1), surface_probes(:,basilisk_probe_number),'DisplayName',titles(i));
    xlabel('t[s]')
    ylabel('surface eleveation[m]')
    xlim([3.2 5.2])
    %mark the period
    %plot(0:T:60,0, 'x')
end
legend();
i = i + 1;

%% plot the piston positions
figure;
hold on
run_number =1;
run_folder = sprintf("%s%d/" , lab_folder, run_number);
load(append(run_folder, "fil3.dat"));
fil = fil3;
plot(t, fil-fil(1,1))


f = 1.425;
xpos=0.042*0.308*max(0,tanh(2*t-0.2).^2).*sin(2*pi*f*t-0.34)*100;


%u = 0.04*0.308*(sech(t).^2.*sin(2*f*pi*t) + 2*f*pi*cos(2*f*pi*t).*tanh(t));
plot(t, xpos);
legend("fil3.dat", "0.042*0.308*max(0,tanh(2*t-0.2).^2).*sin(2*pi*f*t-0.34)*100, f=1.425")

%% plot the piston velocities
figure;
hold on
basilisk_folder = basilisk_folders(3);
surface_probes = readtable(append(basilisk_folder, "surface_probes.csv"), ReadVariableNames=true);
surface_probes = table2array(surface_probes);
plot(surface_probes(:,1), surface_probes(:,2),'DisplayName',"u_x");

plot(t(1:end-1), diff(xpos))
plot(t(1:end-1), diff(fil))
%plot(t, u);
legend(["u_x", "diff(xpos)", "diff(fil3)", "u(t)"])