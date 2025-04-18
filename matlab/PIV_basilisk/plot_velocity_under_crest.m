function [y_scaled, u_crest_scaled] = plot_velocity_under_crest(run_number, pair_number, image_params)
%%plots the velocity under the crest of the wave and returns the velocity
%%and y coordinates under the crest scaled
%% y_scaled = y/water_depth
%% u_crest_scaled = u/(a*omega)

plot_crest_finding = false;


%% load data and parameters
load velocities.mat velocities
load params.mat params

%load parameters for plotting analytical solutions
water_depth = image_params('water_depth');
water_depth = water_depth(run_number)

p = params(run_number);
a = p('a');
a_std = p('std_a');

FentonResult = FentonDispSolver('h', 0.6, 'H', 2*a, 'T', 1/1.425, 'mode', 1);
k = FentonResult.k;

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


%% Find crest

%find the crest by finding the vertical slice with the smallest vertical
%veloctiy
[crest_idx, v_mean, v_mean_lowpass] = find_crest(Vw);
if plot_crest_finding
    figure;
    hold on
    plot(xw(1,:), v_mean_lowpass)
    plot(xw(1,:), v_mean, 'r'); 
    plot(xw(1,crest_idx), v_mean_lowpass(crest_idx), 'xb')
    title('mean v over y')
    legend('lowpass', 'data', 'chosen x point')
end

%% Velocity profile
%%create mask for world based on the three columns closest to the crest 
% coordinates and scale y
% u_crest = Uw(:,crest_idx-1:crest_idx+1);
% u_crest = mean(u_crest, 2);
% u_crest = u_crest(idx(:,crest_idx)&idx(:,crest_idx-1)&idx(:,crest_idx+1));
% crest_mask = idx(:,crest_idx)&idx(:,crest_idx-1)&idx(:,crest_idx+1);
u_crest = Uw(:,crest_idx);
u_crest = u_crest(idx(:,crest_idx)&true);
crest_mask = idx(:,crest_idx)&true;

u_crest_scaled = 1/(a*omega)*u_crest;
y_scaled = 1/water_depth*yw(crest_mask);

%% create shaded area based on the standard deviation of the amplitude'
yw_at_crest = yw(:,crest_idx);

yw_min = min(yw_at_crest):0.0001:a-a_std;
y_min = 1/water_depth*yw_min;

yw_max = min(yw_at_crest):0.0001:a+a_std;
y_max = 1/water_depth*yw_max;


[~, u_max,~,~,~,~,~, ~,~,~,~,~,~] = FentonU(k, water_depth, a+a_std, 0, yw_max);
[~, u_min,~,~,~,~,~, ~,~,~,~,~,~] = FentonU(k, water_depth, a-a_std, 0, yw_min);

u_min = transpose(1/(a*omega)*u_min(:,1));
u_max = transpose(1/(a*omega)*u_max(:,1));

figure;
hold on

f = fill([u_min flip(u_max)], [y_min flip(y_max)], [0.8 0.8 0.8]);
f.EdgeAlpha = 0;
% plot the analytical velocity


yw_analytical = min(yw_at_crest):0.0001:a;
[~, u,~,~,~,~,~, ~,~,~,~,~,~] = FentonU(k, water_depth, a, 0, yw_analytical);
y_analytical_scaled = 1/water_depth*yw_analytical;
plot(1/(a*omega)*u, y_analytical_scaled);


%plot calculated velocity
plot(abs(u_crest_scaled), y_scaled, 'x')


legend('one standard deviation in the amplitude', 'analytical solution', 'experimental','Location','southeast')
title(sprintf('horizontal velocity under the crest run:%d, wave pair:%d', run_number, pair_number))

xlabel('$\frac{v}{a\omega}$', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\frac{y}{h}$', 'interpreter', 'latex', 'FontSize', 20, 'rotation', 0)
%% plot Stokes 5th order
Result = StokesDispSolver('h', 0.6, 'H', 2*a, 'T', 1/1.425, 'mode', 1);
%FentonResult = FentonDispSolver('h', 0.6, 'H', 2*a, 'T', 1/1.425, 'mode', 1);
[~, u,~,~,~,~,~, ~,~,~,~,~,~] = StokesU(Result.k, 0.6, a, 0, yw_analytical);
plot(1/(a*omega)*u, yw_analytical/water_depth, 'DisplayName', 'Zhao')
    
%plot(1/(a*omega)*FentonU(FentonResult.k, 0.6, a*2, 0, yw_analytical), yw_analytical/water_depth, 'DisplayName', 'Fenton')
end
