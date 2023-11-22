%Setup
current_path = cd;
cd '/home/martin/Documents/master/matlab/HydrolabPIV';
setup_hlpiv
cd(current_path)
image_names
coord_config
h = 0.33; %heightof the water surface at rest
n_runs = 15;%number of runs in the wave tank
n_waves = 3;%number of wave image pairs to use per run
close all

%% Perform PIV
force_PIV = false; %true: Force recalculating the PIV or 
%                   false:just PIV on the ones that 
%                         aren't saved in velocities.mat



%%Perform PIV

for run_number=1:n_runs
    for wave=1:n_waves
        disp("run:" + run_number + " | wave:" + wave)
        if force_PIV
            perform_PIV(run_number, wave); %#ok<UNRCH>
        else
            try
                load("velocities.mat")
                run_velocities = velocities(run_number);
                run_velocities(wave);
                disp("PIV already performed moving on to next wave")
            catch
                perform_PIV(run_number, wave);
            end
        end
    end
end

%% Surface height and calculate the amplitude of the waves
for run_number=1:n_runs
    surface_height(run_number);
end
close all

%% Quiver plot
run_number = 8;
pair_number = 1;
max_arrows = 50; %maximum number of arrows in one dimension in the quiver plot
quiver_plot(run_number, pair_number, max_arrows)


%% plot velocity under crest
run_number = 3;
pair_number = 2;
plot_velocity_under_crest(run_number, pair_number);

%% Plot Alpha for one wave 
run_number = 8;
pair_number = 3;
create_plot = true;
plot_alpha(run_number, pair_number, create_plot)

%% Plot mean alpha for the three chosen waves from one run
close all
clear alpa_mean
run_number = 14;
create_plot = false;
for i=1:3
    if i==1
        [alpha_mean, alpha_theoretical, y_at_crest, y_below_crest_scaled] = plot_alpha(run_number , i, create_plot);
    else
        [alpha, ~, ~, ~] = plot_alpha(run_number, i, create_plot);
        alpha_mean = alpha + alpha_mean;
    end
end
alpha_mean = alpha_mean/3;
figure;
hold on
plot(alpha_mean,y_at_crest, 'x')
plot(alpha_theoretical, y_below_crest_scaled)
title(sprintf('mean alpha for run %d', run_number))
xlabel('$\alpha$', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\frac{y}{h}$', 'interpreter', 'latex', 'FontSize', 20, 'rotation', 0)

%% Compare velocity profile to Stokes 5th order
number_of_runs = 15;
load params.mat params

for run_number=1:number_of_runs
    msg = sprintf('calculating 5th order k for run:%d', run_number);
    disp(msg)
    k_5th = find_stokes5th_k(run_number);
    p = params(run_number);
    k = p('k');
    k_str = sprintf("k from second order %f, k from 5th order %f", k, k_5th);
    disp(k_str)
    close all
end


%%
close all
for run_number=1:number_of_runs
    for pair_number=1:n_waves
        compare_stokes_5th_velocity_profile(run_number, pair_number)
    end
end

