function [] = plot_basilisk_alpha(timestep,a, k, omega, h)
%plot alpha = omega/(a*g*k)*(u.^2 + v.^2).^.5
% vs y scaled of the basilisk results
 
%load the basilisk results for the given timestep
timestep_name = sprintf("basilisk_results/timestep_%d", timestep);
load(timestep_name, "U", "X", "mask");
y = X(:,:,2);
u = U(:,:,1);
v = U(:,:,2);

g = 9.81;
alpha = omega/(a*g*k)*(u.^2 + v.^2).^.5;


mean_alpha = zeros(size(alpha, 1), 1)*NaN;
for i=1:size(alpha, 1)
    horizontal_slice_alpha = alpha(i,:);
    mean_alpha(i) = mean(horizontal_slice_alpha(mask(i,:)), 'omitnan');
end


crest_idx = find_crest_index_from_mask(y, mask);
crest_y = y(:,crest_idx);
crest_y = crest_y(mask(:,crest_idx));
plot(mean_alpha(~isnan(mean_alpha)), crest_y/h, 'DisplayName',sprintf('basilisk timestep %d', timestep))
end