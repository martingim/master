% %choose which day to use the results from 
%%
experiment_day = 1;
lab_run_numbers =[1 4 5];
lab_folder = "Lab_results/24_04_12/";
lab_surface_probes = ["f1425_a0_308.csv";
            "f1425_a0_308_r2.csv";
            "f1425_a0_308_r3.csv";
            "f1425_a0_616_r4.csv";  
            "f1425_a0_154_r5.csv"];


%%
experiment_day =2;
lab_run_numbers = [1 4 7]; %runs to plot
sensors = [1 2 4]; %sensors to generate plots for
lab_folder = "Lab_results/24_09_18/";
lab_surface_probes = ["f1425_a308_r1.csv";
            "f1425_a308_r2.csv";
            "f1425_a308_r3.csv";
            "f1425_a616_r4.csv";
            "f1425_a616_r5.csv";
            "f1425_a616_r6.csv";
            "f1425_a154_r7.csv";
            "f1425_a154_r8.csv";
            "f1425_a154_r9.csv"];

%%
%% plot sensor 1 for different runs in different plots
close all
probe_positions = [8 10.05 10.75 11.50;
                   1.5 10.05 10.75 11.50];
%values to start the plots and stokes 5th order fit from
tmin_plot = [20 25 28 30;
        20 25 28 30];


for sensor = sensors
    tmin = tmin_plot(experiment_day, sensor);
    tmax = tmin+5;

    figure('Position', [1000, 818,1000,800]);
    fontsize(25, "points")
    tiledlayout(3,1)
    for lab_run_number=lab_run_numbers
        file = load(append(lab_folder, lab_surface_probes(lab_run_number)));
        t = 0:0.008:(size(file, 1)-1)/125;
        
        nexttile;
        hold on
        %plot the lab results
        col = sensor+2;
        surface_probe = -(file(:,col)-file(1,col));
        max(file(:,col)) - min(file(:,col))/2;
        plot(t, surface_probe, 'DisplayName','surface probe')
        xlabel('t[s]')
        ylabel('surface eleveation[m]')
        xlim([tmin tmax])
        %mark the period
        %plot(0:T:60,0, 'x')
        min_separation = 10;
        lowpass_surface_probe = lowpass(surface_probe, 5, 125);
        surface_probe = surface_probe((t>tmin)&(t<tmax));
        lowpass_surface_probe = lowpass_surface_probe((t>tmin)&(t<tmax));
        
        t_ = t((t>tmin)&(t<tmax));
        LMax = islocalmax(lowpass_surface_probe, 'MinSeparation', min_separation);
        LMin = islocalmin(lowpass_surface_probe, 'MinSeparation', min_separation);
        a = (mean(lowpass_surface_probe(LMax))-mean(lowpass_surface_probe(LMin)))/2;
        %a = (mean(surface_probe(LMax))-mean(surface_probe(LMin)))/2;
        h = 0.6;
        T = 1/1.425;
        Result = StokesDispSolver('h', h, 'H', 2*a, 'T', T, 'mode', 1);
        title(sprintf('Run %d, surface probe at x=%.2fm, ak=%.2f', lab_run_number,probe_positions(experiment_day, sensor), Result.k*a))
        %plot([0 100], [mean(surface_probe(LMax)) mean(surface_probe(LMax))],'--', 'DisplayName', 'mean crest height')
        %plot([0 100], [mean(surface_probe(LMin)) mean(surface_probe(LMin))],'--', 'DisplayName', 'mean trough depth')
        %plot(t, lowpass_surface_probe);
        omega = 2*pi*1.425;
        t_maxima = t_(LMax);
        theta = omega*(t-t_maxima(2));
        [eta, eta1, eta2, eta3,eta4,eta5] = StokesEta(Result.k, h, Result.a, theta);
        plot(t, eta, 'DisplayName','Stokes 5th order analytical, Zhao')
        FentonResult = FentonDispSolver('h', h, 'H', 2*a, 'T', T, 'mode', 1);
        [etaFen, etaFenorder] = FentonEta(FentonResult.k, h, FentonResult.a, theta);
        plot(t, etaFen, 'DisplayName','Stokes 5th order analytical, Fenton')
        %plot(t_(LMax), surface_probe(LMax), 'x', 'DisplayName','Crest location')
        %plot(t_(LMin), surface_probe(LMin), 'x', 'DisplayName','Trough location')
        (2*a)*Result.k;
        legend('Location','southeast')
        [Zhao_crest_height temp]= StokesEta(Result.k, h, Result.a, 0);
        [Fenton_crest_height temp]=FentonEta(FentonResult.k, h, FentonResult.a, 0);
        disp(sprintf("run %d difference in surface height at the crest between Fenont and Zhao:%.2f%%", lab_run_number, (Zhao_crest_height-Fenton_crest_height)/Fenton_crest_height*100));
    end
    print(sprintf('~/Documents/master/movies_and_figures/surface_measured_vs_analytical_sensor%d', sensor), '-dpng')
    
end