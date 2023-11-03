%[SIGMA,T,K,LAMBDA,CP,CG]=WPARAM(F,H)
% This function calculates these different wave parameters using the 
% full dispersion relationship (sigma^2=g*k*tanh(k*H)) given the 
% depth of the water and the f.

function  [sigma,T,k,lambda,Cp,Cg]=wparam(f,H);



g=9.8;
sigma=2*pi*f;
T=1/f;

k0=0.00001;
k1=sigma^2/(g*tanh(k0*H));

while abs(k1-k0)>0.00001
    k0=k1;
    k1=sigma^2/(g*tanh(k0*H));
end

k=k1;
lambda=2*pi/k;
%Cp=sqrt((g*tanh(k*H))/k);
Cp=sigma/k;
Cg=(1/(2*sigma))*(g*tanh(k*H)+g*k*H*(sech(k*H))^2);

