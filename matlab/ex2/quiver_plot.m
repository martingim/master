function [] = quiver_plot(run_number, pair_number, max_arrows)


scale = 3; %for the arrows in the quiver plots
height = 0.33; 

%% load data and parameters
load('velocities.mat')
load('params.mat')

%load parameters for plotting analytical solutions
p = params(run_number);
a = p('a');
k = p('k');
omega = p('omega');
g = 9.82;
t = 0;
%load the velocity fields and world coordinates
run_velocities = velocities(run_number);
UVw = run_velocities(pair_number);
Uw = squeeze(UVw(1,:,:));
Vw = squeeze(UVw(2,:,:));
xw = squeeze(UVw(3,:,:));
yw = squeeze(UVw(4,:,:));
idx = squeeze(UVw(5,:,:));


%%  analytical solution
%potential of wave moving to the right
%phi = a*g/omega*exp(k*y)*sin(k*x-omega*t);
u = @(x, y) a*k*g/omega*exp(k*y).*cos(k*x-omega*t);
v = @(x, y) a*k*g/omega*exp(k*y).*sin(k*x-omega*t);


%% make New mask for quiver plots based on max_arrows
quiver_idx = idx;
if size(quiver_idx, 2)>max_arrows || size(quiver_idx, 1)>max_arrows
    quiver_idx = zeros(size(quiver_idx));
    if size(quiver_idx, 1)>max_arrows
        interval_x = floor(size(quiver_idx, 1)/max_arrows);
    else
        interval_x = 1;
    end
    if size(quiver_idx, 2)>max_arrows
        interval_y = floor(size(quiver_idx, 2)/max_arrows);
    else
        interval_y = 1;
    end
    
    quiver_idx(1:interval_x:end, 1:interval_y:end) = 1;
end
quiver_idx = logical(quiver_idx.*idx);

%% Quiver plots

%World
figure;
hold on;
quiver(xw(quiver_idx),yw(quiver_idx),Uw(quiver_idx),Vw(quiver_idx), scale);
legend('velocity')
title('world')
xlabel('x[m]')
ylabel('y[m]')
hold off;

%Analytical
figure;
quiver(xw(quiver_idx), yw(quiver_idx), u(xw(quiver_idx), yw(quiver_idx)), v(xw(quiver_idx), yw(quiver_idx)), scale)
title('analyctical')
xlabel('x[m]')
ylabel('y[m]')



end