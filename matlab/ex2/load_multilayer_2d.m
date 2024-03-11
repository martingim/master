function [times,results] = load_multilayer_2d(nx, nl)
files = dir(sprintf('basilisk_results/velocities_nx%d_nl%d*.csv', nx, nl));
n_timesteps = numel(files);
times  = zeros(n_timesteps,1);
n_variables = 8;
results = zeros(n_timesteps,n_variables-2, nl, nx);

for n=1:n_timesteps
    filename = sprintf('%s/%s', files(n).folder, files(n).name);
    Tab = readtable(filename);
    array = table2array(Tab);
    times(n) = array(1,1);
    for i=3:n_variables
        results(n,i-2,:,:) = reshape(array(:,i), nl,nx);    
    end
end