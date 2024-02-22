function [] = compare_stokes_5th_velocity_profile(run_number, pair_number, image_params)
%UNTITLED2 Summary of this function goes here
% 
%%
water_depth = image_params('water_depth');
h = water_depth(run_number);
theta=0;
plot_velocity_under_crest(run_number, pair_number, image_params);
load params.mat params
p = params(run_number);
k = p('k5');
a = p('a');
omega = p('omega5');
z = -h:0.01:a;
[phi, u,~,~,~,~,~, w,~,~,~,~,~] = StokesU(k, h, a, theta, z);
plot(u(z>-0.7*h)/a/omega, z(z>-0.7*h)/h, 'b')
legend('PIV', 'theoretical', '5th order Stokes', 'Location','southeast')
end