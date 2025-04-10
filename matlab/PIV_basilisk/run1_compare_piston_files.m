lab_folder = "Lab_results/24_04_12/";
lab_surface_probes = ["f1425_a0_308.csv";
            "f1425_a0_308_r2.csv";
            "f1425_a0_308_r3.csv";
            "f1425_a0_616_r1.csv";  
            "f1425_a0_154_r1.csv"];
run_number = 1;
%% plot
close all;

run_folder = sprintf("%s%d/" , lab_folder, run_number);
load(append(run_folder, "tid.dat"))
t = tid;

figure('Position', [1000, 818,1000,800]);
tiledlayout(2,1)
nexttile;
hold on
title("Piston position files run number 1, 12/04/2024");
load(append(run_folder, "fil1.dat"));
fil = fil1;
%t = 0:0.01:(size(fil1, 1)-1)/100;
plot(t, -(fil-fil(1,1))*0.044, 'DisplayName','fil1.dat*0.042')

load(append(run_folder, "fil2.dat"));
fil = fil2;
plot(t, (fil-fil(1,1))*0.0046, 'DisplayName','fil2.dat*0.0046')

load(append(run_folder, "fil3.dat"));
fil = fil3;
plot(t, (fil-fil(1,1))*0.01, 'k--','DisplayName','fil3.dat*0.01')

load(append(run_folder, "padle_ut.dat"));

plot(t(1:size(padle_ut, 1)), (padle_ut-padle_ut(1))*0.042, 'DisplayName','padle\_ut.dat*0.042')
xlim([0, 3])

plot(t, 0.042*0.308*tanh(t).*sin(2*pi*1.425*(t)), 'DisplayName','0.042\cdot0.308\cdottanh(t)\cdotsin(2\pift)')
legend('Location','southeast')
xlabel('time [s]')
ylabel('piston position [m]')



nexttile;
hold on;
plot(t(1:size(padle_ut, 1)), (padle_ut-padle_ut(1))*0.042, 'DisplayName','padle\_ut.dat')
plot(t, 0.042*0.308*tanh(t).*sin(2*pi*1.425*(t)), '--','DisplayName','tanh(t)*sin(2\pift)')
plot(t(1:4:end)-0.03, (fil(1:4:end)-fil(1,1))*0.01,'xk', 'DisplayName','fil3.dat shifted by 0.03s', 'MarkerSize',10)
legend('Location','southeast')
xlabel('time [s]')
ylabel('piston position [m]')
xlim([10 13])

fontsize(15, "points")
print('~/Documents/master/movies_and_figures/compare_piston_files', '-dpng')
