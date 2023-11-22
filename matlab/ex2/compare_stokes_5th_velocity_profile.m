function [] = compare_stokes_5th_velocity_profile(run_number, pair_number)
%UNTITLED2 Summary of this function goes here
% 
%%
h = 0.33;
theta=0;
plot_velocity_under_crest(run_number, pair_number);
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