function [] = basilisk_moving_piston_velocity_profile(timestep_name,a, omega, h, x_start, x_end)
%plot the velocity under the crest for the moving piston results from
%basilsk

%load the basilisk results for the given timestep
load(timestep_name, "U", "X", "mask");
y = X(:,:,2);
x = X(:,:,1);

mask(X(:,:,1)<x_start)= 0;
mask(X(:,:,1)>x_end)= 0;
crest_idx = find_crest_index_from_mask(y, mask);
x(1,crest_idx)

% figure;
% quiver(X(:,:,1),X(:,:,2),U(:,:,1).*mask,U(:,:,2).*mask)
% figure;
u = U(:,crest_idx,1);
y = X(:,crest_idx,2);
mask = mask(:,crest_idx,1);
omega=8.95;
plot(u(mask)/a/omega, y(mask)/h, 'DisplayName','basilisk')
end