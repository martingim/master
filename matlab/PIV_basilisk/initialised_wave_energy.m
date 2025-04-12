close all;
%% Navier Stokes Solver
close all;
energy_files = [];
legends = [];
energy_files = [energy_files; "~/Documents/master/basilisk/2d_piston/piston-moving/results/LEVEL11_0/energy.csv"]; legends = [legends; "moving piston level 11"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers40/energy_nx128_nl40.csv"]; legends = [legends; "nx:128, l:40s, 377s"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/LEVEL7_nwaves1/energy.txt"]; legends = [legends; "nx 128 dt=0.01"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/LEVEL8_nwaves1/energy.txt"]; legends = [legends; "nx 256 dt=0.005"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/LEVEL9_nwaves1/energy.txt"]; legends = [legends; "nx 512 dt=0.0025"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/dt_small/LEVEL7_nwaves1/energy.txt"]; legends = [legends; "nx 128  dt=0.001"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/dt_small/LEVEL8_nwaves1/energy.txt"]; legends = [legends; "nx 256  dt=0.0005"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/dt_small/LEVEL9_nwaves1/energy.txt"]; legends = [legends; "nx 512  dt=0.00025"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/dt_small/LEVEL10_nwaves1/energy.txt"]; legends = [legends; "nx 1024"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL9_layers30/energy_nx512_nl30.csv"]; legends = [legends; "30 layers, nx 512"];


figure;

hold on;
title("percentage change in Energy compared to k")
xlabel("t [s]")
run = 1;
energy = readtable(energy_files(run));

energy = table2array(energy);
ke0 = energy(1,2);
t = energy(:,1);
ke = energy(:,2)-energy(1,2);
gpe = energy(:,3)-energy(1,3);
te = ke + gpe;
plot(t, ke/ke0*100, 'DisplayName',sprintf('Kinetic Energy'));
plot(t, gpe/ke0*100, 'DisplayName', sprintf('Potential Energy'));
plot(t, (ke+gpe)*100/ke0/2, 'DisplayName', sprintf('Total energy/2'));
title('NS-Solver change in energy');
ylabel('% change in energy');
legend();
%xlim([0,25]);
fontsize(20, "points")
%print('~/Documents/master/movies_and_figures/initialised_NS_energy', '-dpng')


figure;
hold on;
for i =1:size(energy_files)
    energy = readtable(energy_files(i));
    
    energy = table2array(energy);
    ke0 = energy(1,2);
    t = energy(:,1);
    ke = energy(:,2)-energy(1,2);
    gpe = energy(:,3)-energy(1,3);
    te = ke + gpe;
    plot(t, te/ke0*100, 'DisplayName',legends(i))
    
end
xlim([0 10]);
ylabel('% change in energy');
title("Change in  kinetic energy")
legend('Location','southwest')
fontsize(20, "points")
%print('~/Documents/master/movies_and_figures/initialised_NS_energy_LEVEL_comparison', '-dpng')

%% multilayer compare Different LEVELS
close all;

energy_files = [];
legends = [];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers10/energy_nx128_nl10.csv"]; legends = [legends; "10 layers, nx 128"];

energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers8/energy_nx128_nl8.csv"]; legends = [legends; "8 layers, nx 128"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL8_layers15/energy_nx256_nl15.csv"]; legends = [legends; "15 layers, nx 256"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL9_layers30/energy_nx512_nl30.csv"]; legends = [legends; "30 layers, nx 512"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL10_layers60/energy_nx1024_nl60.csv"]; legends = [legends; "60 layers, nx 1024"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers100/energy_nx128_nl100.csv"]; legends = [legends; "100 layers, nx 128"];

figure;

hold on;
title("percentage change in Energy compared to k")
xlabel("t [s]")
run = 1;
energy = readtable(energy_files(run));
energy = table2array(energy);
ke0 = energy(1,2);
t = energy(:,1);
ke = energy(:,2)-energy(1,2);
gpe = energy(:,3)-energy(1,3);
te = ke + gpe;
plot(t, ke/ke0*100, 'DisplayName',sprintf('Kinetic Energy'));
plot(t, gpe/ke0*100, 'DisplayName', sprintf('Potential Energy'));
plot(t, (ke+gpe)*100/ke0/2, 'DisplayName', sprintf('Total energy/2'));
title('Multilayer change in energy');
ylabel('% change in energy');
legend();
xlim([0,25]);
fontsize(20, "points")
print('~/Documents/master/movies_and_figures/initialised_multilayer_energy', '-dpng')


figure;
hold on;
for i =1:size(energy_files)
    energy = readtable(energy_files(i));
    
    energy = table2array(energy);
    ke0 = energy(1,2);
    t = energy(:,1);
    ke = energy(:,2)-energy(1,2);
    gpe = energy(:,3)-energy(1,3);
    te = ke + gpe;
    plot(t, te/ke0*100, 'DisplayName',legends(i))
    
end
xlim([0 25]);
ylabel('% change in energy');
title("Change in  kinetic energy")
legend('Location','southwest')
fontsize(20, "points")
ylim([-12 2]);
print('~/Documents/master/movies_and_figures/initialised_multilayer_energy_LEVEL_comparison', '-dpng')

%% Multilayer compare different layers

close all;

energy_files = [];
legends = [];
%LEVEL7
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers10/energy_nx128_nl10.csv"]; legends = [legends; "nx:128, l:10, 98s"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers20/energy_nx128_nl20.csv"]; legends = [legends; "nx:128, l:20, 216s"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers40/energy_nx128_nl40.csv"]; legends = [legends; "nx:128, l:40s, 377s"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers100/energy_nx128_nl100.csv"]; legends = [legends; "nx:128, l:100, 1210s"];
%LEVEL8
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL8_layers10/energy_nx256_nl10.csv"]; legends = [legends; "nx:256, l:10, 156s"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL8_layers20/energy_nx256_nl20.csv"]; legends = [legends; "nx:256, l:20, 340s"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL8_layers40/energy_nx256_nl40.csv"]; legends = [legends; "nx:256, l:40, 675s"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL8_layers100/energy_nx256_nl100.csv"]; legends = [legends; "nx:256, l:100, 2741s"];
%LEVEL9
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL9_layers20/energy_nx512_nl20.csv"]; legends = [legends; "nx:512, l:20, 626s"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL9_layers40/energy_nx512_nl40.csv"]; legends = [legends; "nx:512, l:40, 1676s"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL9_layers100/energy_nx512_nl100.csv"]; legends = [legends; "nx:512, l:100, 6696s"];
%LEVEL10
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL10_layers40/energy_nx1024_nl40.csv"]; legends = [legends; "nx:1024, l:40, s"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL10_layers10/energy_nx1024_nl10.csv"]; legends = [legends; "nx:1024, l:10, s"];

figure;
hold on;
for i =1:size(energy_files)
    energy = readtable(energy_files(i));
    
    energy = table2array(energy);
    ke0 = energy(1,2);
    t = energy(:,1);
    ke = energy(:,2)-energy(1,2);
    gpe = energy(:,3)-energy(1,3);
    te = ke + gpe;
    plot(t, te/ke0*100, 'DisplayName',legends(i))
    
end
xlim([0 100]);
ylabel('% change in energy');
title("Change in  kinetic energy")
legend('Location','southwest')
fontsize(20, "points")
ylim([-40 2]);
print('~/Documents/master/movies_and_figures/initialised_multilayer_energy_layer_comparison', '-dpng')
