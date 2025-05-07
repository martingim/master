
%% run 5 
%test refineparam
% 
% basilisk_folders = [];
% titles = [];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam1/"];     titles = [titles "LEVEL 14 refineparam 1"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam01/"];    titles = [titles "LEVEL 14 refineparam 0.1"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam001/"];   titles = [titles "LEVEL 14 refineparam 0.01"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam0001/"];  titles = [titles "run 5 LEVEL 14 refineparam 0.001"];

%test levels
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL15/"];  titles = [titles "run 5 LEVEL 15"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam0001/"];  titles = [titles "run 5 LEVEL 14 refineparam 0.001"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL13_refineparam0001/"];  titles = [titles "run 5 LEVEL 13 refineparam 0.001"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL12/"];  titles = [titles "run 5 LEVEL 12"];
%% results for the wave from run number 5 test LEVEL
%half amplitude piston movement
lab_run_number = 5;
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL12/"]; titles = [titles "LEVEL 12"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL13/"]; titles = [titles "LEVEL 13"];
basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14/"]; titles = [titles "LEVEL 14"];
%basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL15/"]; titles = [titles "LEVEL 15"];

%% run 1
basilisk_folders = [];
titles = [];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL11/"];  titles = [titles "run 1 LEVEL 11"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/LEVEL12/"];  titles = [titles "run 1 LEVEL 12"];
basilisk_folders = [basilisk_folders "~/Documents/master/basilisk/2d_piston/boundary-piston/results/run1/conserving/LEVEL12/"];  titles = [titles "conserving run 1 LEVEL 12"];
multilayer_titles = []; multilayer_folders = [];
multilayer_folders = [multilayer_folders "~/Documents/master/basilisk/2d_piston/multilayer-piston/results/run1/LEVEL12_layers20/"]; multilayer_titles = [multilayer_titles "LEVEL12 nl20"];
%%
figure;
hold on;

timestep = 30;
omega = 8.95;
h = 0.6;

mat_files = [];
vtu_files = [];
for i=1:size(basilisk_folders, 2)
    basilisk_folder = basilisk_folders(i);
    mat_files = [mat_files append(basilisk_folder, sprintf("/vtu/matlab/moving_piston_timestep_%d.mat", timestep))];
    vtu_files = [vtu_files append(basilisk_folder, sprintf("/vtu/ascii/%d_0.vtu", timestep))];
    vtu_file = vtu_files(i);
    mat_file = mat_files(i);
    % convert the files from binary vtu to mat files
    if isfile(vtu_file)
        disp("already made ascii file")
    else
        disp("did not find ascii vtu file converting vtu file")
        command = sprintf("/home/martin/ParaView-5.12.0-RC3-MPI-Linux-Python3.10-x86_64/bin/pvpython ~/Documents/master/basilisk/paraview_script.py %s %d", basilisk_folder, timestep);
        re = system(command);
    end
    
    if isfile(mat_file)
        disp("already converted to mat file")
    else
        %convert the vtu files to ascii
        %convert the vtu for the timestep to mat file
        command = sprintf("python3 /home/martin/Documents/master/basilisk/read_vtu_python.py %s %d", basilisk_folder, timestep);
        system(command);
    end
    load(mat_file)
    disp("file loaded")

    % find amplitude using surface_probes.csv
    surface_probes = readtable(append(basilisk_folder, "surface_probes.csv"), 'ReadVariableNames',true, 'VariableNamingRule','preserve');
    probe_names = surface_probes.Properties.VariableNames;
    probe_names{20} = probe_names{20}(1:8);
    basilisk_probe_locations = str2double(probe_names');
    basilisk_probe_locations = basilisk_probe_locations(2:end);
    surface_probes = table2array(surface_probes);
    probe_index = find(abs(basilisk_probe_locations-12.5)<1e-6)+1;

    % surface_probes = load(append(basilisk_folder, "/surface_probes.csv"));
    surface_probes = surface_probes(surface_probes(:,1)>timestep-2,:);
    surface_probes = surface_probes(surface_probes(:,1)<timestep+2,:);
    %plot(surface_probes(:,2))
    min_separation = 10;
    LMax = islocalmax(surface_probes(:,probe_index), 'MinSeparation', min_separation);
    LMin = islocalmin(surface_probes(:,probe_index), 'MinSeparation', min_separation);
    a = (mean(surface_probes(LMax,probe_index))-mean(surface_probes(LMin,probe_index)))/2
    
    % plot the velocity profile
    x_start = 12.1; %where to plot the mean kinetic energy from and to
    x_end = 12.9;
    
    
    mask(X(:,:,1)<x_start) = 0;
    mask(X(:,:,1)>x_end) = 0;
    y = X(:,:,2);
    u = U(:,:,1);
    v = U(:,:,2);
    [temp, m_idx] = max(sum(mask, 1));
    X(1,m_idx,1)
    % figure
    % plot(sum(mask, 1))
    
    for m_idx =m_idx:m_idx
        plot(u(mask(:,m_idx),m_idx), y(mask(:,m_idx),m_idx)/h, "DisplayName",titles(i))
    end
    z = -h:0.001:a;
    T =2*pi/omega;
    Result = StokesDispSolver('h', h, 'H', 2*a, 'T', T, 'mode', 1);
    a
    Result.a
    [~, stku,~,~,~,~,~, ~,~,~,~,~,~] = StokesU(Result.k, h, a, 0, z);
    plot(stku, z/h, 'DisplayName',"Zhao")
    legend
    % figure
    % 
    % g = 9.81;
    % alpha = omega/(a*g*k)*(u.^2 + v.^2).^.5;
    % mean_alpha = zeros(size(alpha, 1), 1)*NaN;
    % for j=1:size(alpha, 1)
    %     horizontal_slice_alpha = alpha(j,:);
    %     mean_alpha(j) = mean(horizontal_slice_alpha(mask(j,:)), 'omitnan');
    % end
    % crest_idx = find_crest_index_from_mask(y, mask);
    % 
    % mean_alpha = mean_alpha(mask(:,crest_idx));
    % crest_y = y(:,crest_idx);
    % crest_y = crest_y(mask(:,crest_idx));
    % plot(mean_alpha, crest_y/h, 'DisplayName',titles(i))
    % legend()
    % title("horizontal mean of velocity")
    % xlabel('$\alpha=\frac{\omega}{agk}(u^2+v^2)^{\frac{1}{2}}$', 'interpreter', 'latex', 'FontSize', 20)
    % ylabel('$\frac{y}{h}$', 'interpreter', 'latex', 'FontSize', 20, 'rotation', 0)
end
%%
multilayer_mat_files = [];

for i=1:size(multilayer_folders, 2)
    basilisk_folder = multilayer_folders(i);
    multilayer_mat_files = [multilayer_mat_files append(basilisk_folder, sprintf("/vts/matlab/timestep_%d.mat", timestep))];
    mat_file = multilayer_mat_files(i);
    %% convert the files from binary vtu to mat files
    
    if isfile(mat_file)
        disp("already converted to mat file")
    else
        %convert the vtu files to ascii
        %convert the vtu for the timestep to mat file
        command = sprintf("python3 /home/martin/Documents/master/basilisk/read_vts_python.py %s %d", basilisk_folder, timestep);
        system(command);
    end
    load(mat_file)
    disp("file loaded")

    % find amplitude using surface_probes.csv
    surface_probes = readtable(append(basilisk_folder, "surface_probes.csv"), 'ReadVariableNames',true, 'VariableNamingRule','preserve');
    probe_names = surface_probes.Properties.VariableNames;
    probe_names{20} = probe_names{20}(1:8);
    basilisk_probe_locations = str2double(probe_names');
    basilisk_probe_locations = basilisk_probe_locations(2:end);
    surface_probes = table2array(surface_probes);
    probe_index = find(abs(basilisk_probe_locations-12.5)<1e-6)+1;

    % surface_probes = load(append(basilisk_folder, "/surface_probes.csv"));
    surface_probes = surface_probes(surface_probes(:,1)>timestep-2,:);
    surface_probes = surface_probes(surface_probes(:,1)<timestep+2,:);
    %plot(surface_probes(:,2))
    min_separation = 10;
    LMax = islocalmax(surface_probes(:,probe_index), 'MinSeparation', min_separation);
    LMin = islocalmin(surface_probes(:,probe_index), 'MinSeparation', min_separation);
    a = (mean(surface_probes(LMax,probe_index))-mean(surface_probes(LMin,probe_index)))/2
    
    % plot the velocity profile
    x_start = 12; %where to plot the mean kinetic energy from and to
    x_end = 13;
    
    x = X(:,:,1);
    y = X(:,:,2);
    y(x<x_start) = 0;
    y(x>x_end) = 0;
    u = U(:,:,1);
    v = U(:,:,2);
    [temp, m_idx] = max(y(1,:));
    
    % figure
    % plot(y(1,:))
    
    for m_idx =m_idx
        plot(u(:,m_idx), y(:,m_idx)/h, "DisplayName",multilayer_titles(i))
    end
    z = -h:0.001:a;
    T =2*pi/omega;
    Result = StokesDispSolver('h', h, 'H', 2*a, 'T', T, 'mode', 1);
    a
    Result.a
    [~, stku,~,~,~,~,~, ~,~,~,~,~,~] = StokesU(Result.k, h, a, 0, z);
    plot(stku, z/h, 'DisplayName',"Zhao")
    legend
    % figure
    % 
    % g = 9.81;
    % alpha = omega/(a*g*k)*(u.^2 + v.^2).^.5;
    % mean_alpha = zeros(size(alpha, 1), 1)*NaN;
    % for j=1:size(alpha, 1)
    %     horizontal_slice_alpha = alpha(j,:);
    %     mean_alpha(j) = mean(horizontal_slice_alpha(mask(j,:)), 'omitnan');
    % end
    % crest_idx = find_crest_index_from_mask(y, mask);
    % 
    % mean_alpha = mean_alpha(mask(:,crest_idx));
    % crest_y = y(:,crest_idx);
    % crest_y = crest_y(mask(:,crest_idx));
    % plot(mean_alpha, crest_y/h, 'DisplayName',titles(i))
    % legend()
    % title("horizontal mean of velocity")
    % xlabel('$\alpha=\frac{\omega}{agk}(u^2+v^2)^{\frac{1}{2}}$', 'interpreter', 'latex', 'FontSize', 20)
    % ylabel('$\frac{y}{h}$', 'interpreter', 'latex', 'FontSize', 20, 'rotation', 0)
end
