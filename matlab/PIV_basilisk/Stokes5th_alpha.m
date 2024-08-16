function [k, alpha, alpha_3, z] = Stokes5th_alpha(a, omega, h, plot_figure)
%a = measured amplitude of the wave
% omega = frequency of the wave
% h = the water level in the tank
g = 9.81;
H = a*2;
T = 2*pi/omega;


Result = StokesDispSolver('h', h, 'H', H,  'T', T, 'mode', 1);


z = linspace(-h,a, 1000);
[phi, u,u1,u2,u3,~,~, w,w1,w2,w3,~,~] = StokesU(Result.k, Result.h, Result.a, 0, z);
u_3 = u1+u2+u3;
w_3 = w1+w2+w3; 
alpha = omega/(a*g*Result.k)*(u.^2+w.^2).^.5;
alpha_3 = omega/(a*g*Result.k)*(u_3.^2+w_3.^2).^.5;
alpha_1 = omega/(a*g*Result.k)*(u1.^2+w1.^2).^.5;
if plot_figure
    hold on
    plot(alpha, z/h, 'DisplayName',"Stokes 5th theoretical")
    plot(alpha_3, z/h, 'DisplayName',"Stokes 3rd theoretical")
    plot(alpha_1, z/h, 'DisplayName', "Stokes linear theoretical")
end
k = Result.k;
end