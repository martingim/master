function [] = basilisk_moving_piston_velocity_profile(timestep,a, omega, h)
%plot the velocity under the crest for the moving piston results from
%basilsk

%load the basilisk results for the given timestep
timestep_name = sprintf("basilisk_results/moving_piston_timestep_%d", timestep);
load(timestep_name, "U", "X", "mask");
y = X(:,:,2);

mask(:,1:20)= 0;
crest_idx = find_crest_index_from_mask(y, mask);

% figure;
% quiver(X(:,:,1),X(:,:,2),U(:,:,1).*mask,U(:,:,2).*mask)
% figure;
u = U(:,crest_idx,1);
y = X(:,crest_idx,2)-0.6;
mask = mask(:,crest_idx,1);

plot(u(mask)/a/omega, y(mask)/h, 'DisplayName',sprintf('basilisk timestep %d', timestep))
end