clear all
close all
figure;
hold on;
load Lab_results/24_04_12/1/fil1.dat
fil = fil1;
t = 0.008:0.01:size(fil1, 1)/100;

plot(t, fil-fil(1,1))
load Lab_results/24_04_12/1/fil2.dat
fil = fil2;
plot(t, fil-fil(1,1))

load Lab_results/24_04_12/1/fil3.dat
fil = fil3;
plot(t, fil-fil(1,1))

load Lab_results/24_04_12/1/padle_ut.dat
plot(t, (padle_ut-padle_ut(1)))

legend('fil1', 'fil2', 'fil3', 'padle ut')

hs = 0:1/1.425:60;
xs = zeros(size(hs));
plot(hs,xs, 'X');

%%
close all
figure;
hold on
load Lab_results/24_04_12/f1425_a0_308.csv
file = f1425_a0_308;
t = 0.008:0.008:size(file, 1)/125;
col = 7;
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

load basilisk_results/surface_probes.csv
plot(surface_probes(:,1), surface_probes(:,2),'DisplayName','NS')

 
% hs = 0:1/1.425:60;
% xs = zeros(size(hs));
% plot(hs,xs, 'X');
%%
load Lab_results/24_04_12/1/padle_ut.dat
load Lab_results/24_04_12/1/fil3.dat
plot((padle_ut-padle_ut(1)))
