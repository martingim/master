function [] = basilisk_multilayer_alpha(a, k, omega, h)

[times,res] = load_multilayer_2d('basilisk_results/velocities.csv');
times = [times(1) times(end)]; %only plot first and last timestep
g = 9.81;
N = 1000;
for i = 1:size(times,1)
    t = times(i); 
    X = squeeze(res(i,1,:,:));
    Y = squeeze(res(i,2,:,:));
    U = squeeze(res(i,3,:,:));
    V = squeeze(res(i,4,:,:));
    eta=squeeze(res(i,5,:,:));
    alpha = omega/(a*g*k)*(U.^2+V.^2).^.5;
    % calculate the mean
    [ymax, ind] = max(max(Y));
    y = linspace(ymax, -h, N);
    mean_alpha = zeros(N, size(alpha, 2))*NaN;
    for i = 1:size(alpha, 2)
        height = floor((max(Y(:,i))+h)/(ymax+h)*N);
        mean_alpha(N-height+1:end, i) = interp1(Y(:,i),alpha(:,i),y(N-height+1:end));
    end
    mean_alpha = mean(mean_alpha, 2, 'omitnan');
    plot(mean_alpha, y/h, 'DisplayName', sprintf("multilayer at time t=%.1f", t))
end

