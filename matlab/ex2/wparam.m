function  [sigma,T,k,lambda,Cp,Cg]=wparam(f,H)
%[SIGMA,T,K,LAMBDA,CP,CG]=WPARAM(F,H)
% This function calculates different wave parameters using the 
% Second order dispersion relation (sigma=omega.^2/(g*k) and sigma = tanh(k*H)) given the 
% depth of the water H and the frequency f.

g=9.81;
T=1/f;
omega = 2*pi*f;

%% solve for k

syms K;

k = vpasolve(omega.^2/(g*K)-tanh(K*H)==0, K, 5);
%% calculate remaining parameters

sigma = tanh(k*H);
lambda=2*pi/k;

Cp=sigma/k;
Cg=(1/(2*sigma))*(g*tanh(k*H)+g*k*H*(sech(k*H))^2);