function [] = plot_basilisk_velocity_profile(timestep, a, omega, h)
%
 
%load the basilisk results for the given timestep
timestep_name = sprintf("basilisk_results/timestep_%d", timestep);
load(timestep_name, "U", "X", "mask");
y = X(:,:,2);


crest_idx = find_crest_index_from_mask(y, mask);

u = U(:,crest_idx,1);
y = X(:,crest_idx,2);
mask = mask(:,crest_idx,1);
plot(u(mask)/a/omega, y(mask)/h, 'DisplayName',sprintf('basilisk timestep %d', timestep))

end