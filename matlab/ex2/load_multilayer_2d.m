function [times,results] = load_multilayer_2d(filename)

Tab = readtable(filename);
array = table2array(Tab);
nl = max(array(:,2)) + 1; %number of layers
times = unique(array(:,1));
n_timesteps = size(times,1);
n_variables = size(array,2);
nx = numel(array)/n_variables/n_timesteps/nl;
results = zeros(n_timesteps,n_variables-2, nl, nx);
for j=1:n_timesteps
    results_ = array(array(:,1,1)==times(j),3:end);
    for i=1:n_variables-2
        results(j,i,:,:) = reshape(results_(:,i), nl,nx);    
    end
end

end