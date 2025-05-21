close all;
%% Navier Stokes Solver compare energy intialised waves periodic boundary
close all;
energy_files = [];
legends = [];
% energy_files = [energy_files; "~/Documents/master/basilisk/2d_piston/piston-moving/results/LEVEL11_0/energy.csv"]; legends = [legends; "moving piston level 11"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL9_layers40/energy_nx512_nl40.csv"]; legends = [legends; "multilayer nx:512, nl:40s, s"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/CFL0_5/LEVEL7_nwaves1/energy.txt"]; legends = [legends; "CFL:0.5, nx:128, runtime:14s"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/CFL0_5/LEVEL8_nwaves1/energy.txt"]; legends = [legends; "CFL:0.5, nx 256 runtime:43"];

% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/LEVEL7_nwaves1/energy.txt"]; legends = [legends; "nx 128 runtime:14s"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/LEVEL8_nwaves1/energy.txt"]; legends = [legends; "nx 256 runtime:43"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/LEVEL9_nwaves1_dt0015/energy.txt"]; legends = [legends; "maxDT~dx, adaptive, nx:512 runtime:232s"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/LEVEL9_nwaves1_dt0031/energy.txt"]; legends = [legends; "maxDT~2dx, adaptive, nx:512 runtime:143s"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS/results/LEVEL7/energy.txt"]; legends = [legends; "nx:128, runtime:21s"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS/results/LEVEL8/energy.txt"]; legends = [legends; "nx 256, runtime:239"];

%rerun with same dt
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS/results/LEVEL9_dt0015/energy.txt"]; legends = [legends; "maxDT~dx, nx 512, runtime:3264"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS/results/LEVEL9/energy.txt"]; legends = [legends; "maxDT~2dx, nx 512, runtime:2168"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/dt_small/LEVEL7_nwaves1/energy.txt"]; legends = [legends; "nx 128  dt=0.001"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/dt_small/LEVEL8_nwaves1/energy.txt"]; legends = [legends; "nx 256  dt=0.0005"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/dt_small/LEVEL9_nwaves1/energy.txt"]; legends = [legends; "nx 512  dt=0.00025"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/NS-adaptive/results/dt_small/LEVEL10_nwaves1/energy.txt"]; legends = [legends; "nx 1024"];
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
xlabel("t [s]")
legend();
%xlim([0,25]);
fontsize(20, "points")
%print('~/Documents/master/movies_and_figures/initialised_NS_energy', '-dpng')


figure;
hold on;
for i =1:size(energy_files,1)
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
xlabel("t [s]")

title("Change in  kinetic energy")
legend('Location','northwest')
fontsize(20, "points")
print('~/Documents/master/movies_and_figures/initialised_NS_adaptive_multilayer_comp', '-dpng')
%print('~/Documents/master/movies_and_figures/initialised_NS_energy_LEVEL_comparison', '-dpng')

%% multilayer compare Different LEVELS
close all;

energy_files = [];
legends = [];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers10/energy_nx128_nl10.csv"]; legends = [legends; "10 layers, nx 128"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers20/energy_nx128_nl20.csv"]; legends = [legends; "20 layers, nx 128"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers40/energy_nx128_nl40.csv"]; legends = [legends; "40 layers, nx 128"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL8_layers10/energy_nx256_nl10.csv"]; legends = [legends; "10 layers, nx 256"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL8_layers20/energy_nx256_nl20.csv"]; legends = [legends; "20 layers, nx 256"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL8_layers40/energy_nx256_nl40.csv"]; legends = [legends; "40 layers, nx 256"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL9_layers10/energy_nx512_nl10.csv"]; legends = [legends; "10 layers, nx 512"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL9_layers20/energy_nx512_nl20.csv"]; legends = [legends; "20 layers, nx 512"];
energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL9_layers40/energy_nx512_nl40.csv"]; legends = [legends; "40 layers, nx 512"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL10_layers60/energy_nx1024_nl60.csv"]; legends = [legends; "60 layers, nx 1024"];
% energy_files = [energy_files; "~/Documents/master/basilisk/initialised_wave/multilayer/results/LEVEL7_layers100/energy_nx128_nl100.csv"]; legends = [legends; "100 layers, nx 128"];

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
xlim([0 10]);
ylabel('% change in energy');
title("Change in  kinetic energy")
legend('Location','southwest')
fontsize(20, "points")
ylim([-6 2]);
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



%% compare short wave signal energy in large wavetank
close all
energy_files = [];
legends = [];
% 24.6/2^13;

energy_files = [energy_files; "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run11/LEVEL11/energy.csv"]; legends = [legends; "NS boundary dx:12.0mm, runtime:1002"];
energy_files = [energy_files; "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run11/LEVEL12_noDT/energy.csv"]; legends = [legends; "NS boundary dx:6.0mm, runtime:2355"];
energy_files = [energy_files; "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run11/LEVEL12/energy.csv"]; legends = [legends; "NS boundary max dt:0.0025s, dx:6.0mm, runtime:5143"];

energy_files = [energy_files; "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run11/LEVEL12_layers20_noDT/energy.csv"]; legends = [legends; "multilayer dx:6.0mm, nl:20, 805s"];
% energy_files = [energy_files; "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run11/LEVEL12_layers20/energy.csv"]; legends = [legends; "multilayer dx:6.0mm, nl:20, 805s"];
energy_files = [energy_files; "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run11/LEVEL13_layers20/energy.csv"]; legends = [legends; "multilayer dx:3.0mm, nl:20, 2331s"];
                                                 
% energy_files = [energy_files; "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run11/LEVEL12/energy.csv"]; legends = [legends; "boundary dx:, nl:40s, s"];

energy_files = [energy_files; "~/Documents/master/basilisk/2d_piston/piston-moving/results/run11/LEVEL11_0/energy.csv"]; legends = [legends; "moving piston dt~dx:12.0mm"];



comparison_time = 3.5;
comparison_idx = round(comparison_time*100);
figure;

hold on;
title("percentage change in Energy compared to k")
xlabel("t [s]")
run = 1;
energy = readtable(energy_files(run));

energy = table2array(energy);
ke0 = energy(comparison_idx,2);
gpe0 = energy(comparison_idx,3);
t = energy(:,1)-comparison_time;
ke = energy(:,2)-ke0;
gpe = energy(:,3)-gpe0;
te = ke + gpe;
plot(t, ke/ke0*100, 'DisplayName',sprintf('Kinetic Energy'));
plot(t, gpe/ke0*100, 'DisplayName', sprintf('Potential Energy'));
plot(t, (ke+gpe)*100/ke0/2, 'DisplayName', sprintf('Total energy/2'));
title('NS-Solver change in energy');
ylabel('% change in energy');
legend();
xlim([0,25-comparison_time]);
fontsize(20, "points")
%print('~/Documents/master/movies_and_figures/initialised_NS_energy', '-dpng')


figure;
hold on;
for i =1:size(energy_files, 1)
    energy = readtable(energy_files(i));
    
    energy = table2array(energy);
    ke0 = energy(comparison_idx, 2);
    gpe0 = energy(comparison_idx, 3);
    t = energy(:,1);
    ke = energy(:,2)-ke0;
    gpe = energy(:,3)-gpe0;
    te = ke + gpe;
    plot(t, ke/ke0*100, 'DisplayName',legends(i))
    
end
xlim([0 20]);
ylabel('% change in energy');
title("Change in  kinetic energy")
legend('Location','northwest')
fontsize(20, "points")
% print('~/Documents/master/movies_and_figures/initialised_NS_adaptive_multilayer_comp', '-dpng')
% print('~/Documents/master/movies_and_figures/initialised_NS_energy_LEVEL_comparison', '-dpng')