function [y_scaled, u_crest_scaled] = plot_velocity_under_crest(run_number, pair_number)
%%plots the velocity under the crest of the wave and returns the velocity
%%and y coordinates under the crest scaled
%% y_scaled = y/height
%% u_crest_scaled = u/(a*omega)

height = 0.33; 

%% load data and parameters
load velocities.mat velocities
load params.mat params

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


%% Velocity profile

%find the crest by finding the vertical slice with the smallest vertical
%veloctiy
figure;
V = Vw;
V(isnan(Uw))=0;
v_abs = abs(V);
v_mean = mean(v_abs, 1);

%perform lowpass to smooth the data
v_mean_lowpass = lowpass(v_mean, 0.001);
plot(xw(1,:), v_mean_lowpass)
%remove the leftmost and rightmost part of the data
v_mean_lowpass(1:floor(size(v_mean_lowpass,2)*0.25)) = NaN;
v_mean_lowpass(floor(size(v_mean_lowpass,2)*0.75):end) = NaN;

hold on
[ min_v, crest_idx] = min(v_mean_lowpass);
plot(xw(1,:), v_mean, 'r'); 
plot(xw(1,crest_idx), min_v, 'xb')
title('mean v over y')
legend('lowpass', 'data', 'chosen x point')
%find the velocity profile at the crest

u_crest = Uw(:,crest_idx-1:crest_idx+1);
u_crest = mean(u_crest, 2);
u_crest = u_crest(idx(:,crest_idx)&idx(:,crest_idx-1)&idx(:,crest_idx+1));
crest_mask = idx(:,crest_idx)&idx(:,crest_idx-1)&idx(:,crest_idx+1);
u_crest_scaled = 1/(a*omega)*u_crest;
y_scaled = 1/height*yw(crest_mask);

figure;
hold on 
plot(u_crest_scaled, y_scaled, 'x')
plot(1/(a*omega)*u(xw(crest_mask)*0, yw(crest_mask)), y_scaled);

legend('measured', 'theoretical')
title('theoretical and measured horizontal velocity under the crest')

xlabel('$\frac{v}{a\omega}$', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\frac{y}{h}$', 'interpreter', 'latex', 'FontSize', 20, 'rotation', 0)

end
