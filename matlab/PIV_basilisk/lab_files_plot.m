clear all
figure;
hold on;
load Lab_results/24_04_12/1/fil1.dat
fil = fil1;
plot(fil-fil(1,1))

load Lab_results/24_04_12/1/fil2.dat
fil = fil2;
plot(fil-fil(1,1))

load Lab_results/24_04_12/1/fil3.dat
fil = fil3;
plot(fil-fil(1,1))

load Lab_results/24_04_12/1/padle_ut.dat
plot((padle_ut-padle_ut(1)))

legend('fil1', 'fil2', 'fil3', 'padle ut')


%%
figure;
hold on
load Lab_results/24_04_12/f1425_a0_308.csv
file = f1425_a0_308;
t = 0.01:0.01:size(file, 1)/100;
col = 7;
for col=1:6
    max(file(:,col)) - min(file(:,col))/2

    plot(t, file(:,col)-file(1,col))
end

legend('sensor1', 'sensor2', 'sensor3', 'sensor4', 'sensor5', 'sensor6')

load basilisk_results/surface_probes.csv
plot(surface_probes(:,1), surface_probes(:,2))
%%
load Lab_results/24_04_12/1/padle_ut.dat
load Lab_results/24_04_12/1/fil3.dat
plot((padle_ut-padle_ut(1)))
