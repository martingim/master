% This is a solver for nonlinear dispersion relationship of Stokes theory.
% The full nonlinear dispersion relationship has h (mean water depth), either T
% (wave period) or omega0 (leading order frequency), and 
% either a or H is given. Results include k, H, a, omega, omega0. 
% Mode 1, the first order wave amplitude/height is given as a, H.
% Mode 2, the first harmonic term's coefficient is given as aw, Hw. 
% Example 1, 
% Results =StokesDispSolver('h',30,'T', 5, 'a',2.5, 'mode',1)
% Example 2, 
% Results =StokesDispSolver('h',30,'T', 5, 'a',2.5, 'mode',2)
% Example 3, 
% Results =StokesDispSolver('h',30,'T', 5, 'H',5, 'mode',1)
% Example 4, 
% Results =StokesDispSolver('h',1,'omega0', 1.2566, 'a',0.15, 'mode',1)
% Example 5, 
% Results =StokesDispSolver('h',1,'omega0', 1.2566, 'H',0.3, 'mode',1)

function [Results] = StokesDispSolver(varargin)
names = varargin(1:2:end);
values = varargin(2:2:end);
for kIn = 1:numel(names)
    switch names{kIn}
        case 'a'
            a = values{kIn};
        case 'H'
            H = values{kIn};
        case 'h'
            h = values{kIn};
        case 'T'
            T = values{kIn};
        case 'omega0'
            omega0 = values{kIn};
        case 'mode'
            modeNo = values{kIn};
    end
end
Results.h = h;
syms k
g = 9.81; 
alpha1 = cosh(2*k*h);
sigma = tanh(k*h);
if and(exist('a', 'var') ==1, exist('T', 'var') ==1)
    omega0 = sqrt(g*k*tanh(k*h));
    Results.T = T;
    omega = 2*pi/T;
    Results.omega = omega;
    L0 = g*T^2/2/pi;
    if modeNo == 1 
        omega2 = k^2*a^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
        omega4 = (a^4*k^4*(20*alpha1^5 + 112*alpha1^4 ...
        - 100*alpha1^3 - 68*alpha1^2 - 211*alpha1 + 328))/(32*(alpha1 - 1)^5);
        omegaFun = omega0*(1+omega2+omega4) - omega;  
        kResults = double(vpasolve(omegaFun, k, [0 Inf]));
        Results.a = a;
        sigma = tanh(kResults*h);
        B31 = (3+8*sigma^2-9*sigma^4)/16/sigma^4;
        B51 = (121*cosh(2*h*kResults)^5 + 263*cosh(2*h*kResults)^4 +376*cosh(2*h*kResults)^3 ...
        - 1999*cosh(2*h*kResults)^2 + 2509*cosh(2*h*kResults) - 1108)/(192*(cosh(2*h*kResults) - 1)^5);
        Results.aw = a*(1+kResults^2.*a^2*B31+kResults^4.*a^4*B51);
    elseif modeNo == 2 % in case first harmonic wave amplitude aw is given
        syms a0

        B31 = (3+8*sigma^2-9*sigma^4)/16/sigma^4;
        B51 = (121*cosh(2*h*k)^5 + 263*cosh(2*h*k)^4 +376*cosh(2*h*k)^3 ...
        - 1991*cosh(2*h*k)^2 + 2509*cosh(2*h*k) - 1108)/(192*(cosh(2*h*k) - 1)^5);
        a0Fun = a - a0*(1+k^2.*a0^2*B31+k^4.*a0^4*B51);
        omega2 = k^2*a0^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
        omega4 = (a0^4*k^4*(20*alpha1^5 + 112*alpha1^4 ...
        - 100*alpha1^3 - 68*alpha1^2 - 211*alpha1 + 328))/(32*(alpha1 - 1)^5);
        omegaFun = omega0*(1+omega2+omega4) - omega; 
        S = vpasolve([a0Fun , omegaFun], [k, a0],[2*pi/L0,a]);
        kResults = abs(double(S.k));
        Results.a = abs(double(S.a0));
        Results.aw = a;
        a = Results.a;
    end
    Results.k = kResults;
    Results.omega0 = sqrt(g*kResults*tanh(kResults*h));
    sigma = tanh(kResults*h);
    B31 = (3+8*sigma^2-9*sigma^4)/16/sigma^4;
    B33 = (27-9*sigma^2+9*sigma^4-3*sigma^6)/64/sigma^6;   
    B51 = (121*cosh(2*h*kResults)^5 + 263*cosh(2*h*kResults)^4 +376*cosh(2*h*kResults)^3 ...
        - 1999*cosh(2*h*kResults)^2 + 2509*cosh(2*h*kResults) - 1108)/(192*(cosh(2*h*kResults) - 1)^5);
    B53 = (57*cosh(2*h*kResults)^7 + 204*cosh(2*h*kResults)^6 - 53*cosh(2*h*kResults)^5 ...
        -782*cosh(2*h*kResults)^4 - 741*cosh(2*h*kResults)^3 - 52*cosh(2*h*kResults)^2 ...
        + 371*cosh(2*h*kResults) + 186)*9/((3*cosh(2*h*kResults) + 2)*(cosh(2*h*kResults) - 1)^6*128);
    B55 = (300*cosh(2*h*kResults)^8 + 1579*cosh(2*h*kResults)^7 + 3176*cosh(2*h*kResults)^6 + 2949*cosh(2*h*kResults)^5 ...
        + 1188*cosh(2*h*kResults)^4 + 675*cosh(2*h*kResults)^3 + 1326*cosh(2*h*kResults)^2 + 827*cosh(2*h*kResults) ...
        + 130)*5/((cosh(2*h*kResults) - 1)^6*(12*cosh(2*h*kResults)^2 + 11*cosh(2*h*kResults) + 2)*384);
    Hout = 2*a+2*(B31+B33)*kResults^2*a^3+2*(B51+B53+B55)*kResults^4*a^5;
    Results.H = Hout;
 
elseif and(exist('H', 'var') ==1, exist('T', 'var') ==1)
    omega0 = sqrt(g*k*tanh(k*h));
    Results.T = T;
    omega = 2*pi/T;
    Results.omega = omega;
    L0 = g*T^2/2/pi;
    syms a
     sigma = tanh(k*h);
        B33 = (27-9*sigma^2+9*sigma^4-3*sigma^6)/64/sigma^6;        
        B55 = (300*cosh(2*h*k)^8 + 1579*cosh(2*h*k)^7 + 3176*cosh(2*h*k)^6 + 2949*cosh(2*h*k)^5 ...
            + 1188*cosh(2*h*k)^4 + 675*cosh(2*h*k)^3 + 1326*cosh(2*h*k)^2 + 827*cosh(2*h*k) ...
            + 130)*5/((cosh(2*h*k) - 1)^6*(12*cosh(2*h*k)^2 + 11*cosh(2*h*k) + 2)*384);
        B31 = (3+8*sigma^2-9*sigma^4)/16/sigma^4;
        B51 = (121*cosh(2*h*k)^5 + 263*cosh(2*h*k)^4 +376*cosh(2*h*k)^3 ...
        - 1999*cosh(2*h*k)^2 + 2509*cosh(2*h*k) - 1108)/(192*(cosh(2*h*k) - 1)^5);
        B53 = (57*cosh(2*h*k)^7 + 204*cosh(2*h*k)^6 - 53*cosh(2*h*k)^5 ...
            -782*cosh(2*h*k)^4 - 741*cosh(2*h*k)^3 - 52*cosh(2*h*k)^2 ...
            + 371*cosh(2*h*k) + 186)*9/((3*cosh(2*h*k) + 2)*(cosh(2*h*k) - 1)^6*128);
        Results.H = H;
        Hfun1 = H -(2*a+2*(B31+B33)*k^2*a^3+2*(B51+B53+B55)*k^4*a^5) ==0;
        omega2 = k^2*a^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
        omega4 = (a^4*k^4*(20*alpha1^5 + 112*alpha1^4 ...
        - 100*alpha1^3 - 68*alpha1^2 - 211*alpha1 + 328))/(32*(alpha1 - 1)^5);
        omegaFun1 = omega0*(1+omega2+omega4) - omega ==0;
        S = vpasolve([Hfun1 , omegaFun1], [k, a],[2*pi/L0,H/2]);
        kResults = abs(double(S.k));
        aResults = double(S.a);
        Results.a = aResults;
        sigma = tanh(kResults*h);
        B31 = (3+8*sigma^2-9*sigma^4)/16/sigma^4;
        B51 = (121*cosh(2*h*kResults)^5 + 263*cosh(2*h*kResults)^4 +376*cosh(2*h*kResults)^3 ...
        - 1999*cosh(2*h*kResults)^2 + 2509*cosh(2*h*kResults) - 1108)/(192*(cosh(2*h*kResults) - 1)^5);
        Results.aw = aResults*(1+kResults^2.*aResults^2*B31+kResults^4.*aResults^4*B51);
        Results.omega0 = sqrt(g*kResults*tanh(kResults*h));
elseif  and(exist('a', 'var') ==1, exist('omega0', 'var') ==1)
    Results.omega0 = omega0;
    kResults = double(vpasolve(omega0^2 - g*k*tanh(k*h), k, [0 inf]));
    alpha1 = cosh(2*kResults*h);   
    sigma = tanh(kResults*h);
    B33 = (27-9*sigma^2+9*sigma^4-3*sigma^6)/64/sigma^6;    
    B55 = (300*cosh(2*h*kResults)^8 + 1579*cosh(2*h*kResults)^7 + 3176*cosh(2*h*kResults)^6 + 2949*cosh(2*h*kResults)^5 ...
        + 1188*cosh(2*h*kResults)^4 + 675*cosh(2*h*kResults)^3 + 1326*cosh(2*h*kResults)^2 + 827*cosh(2*h*kResults) ...
        + 130)*5/((cosh(2*h*kResults) - 1)^6*(12*cosh(2*h*kResults)^2 + 11*cosh(2*h*kResults) + 2)*384);
    B31 = (3+8*sigma^2-9*sigma^4)/16/sigma^4;
    B51 = (121*cosh(2*h*kResults)^5 + 263*cosh(2*h*kResults)^4 +376*cosh(2*h*kResults)^3 ...
        - 1999*cosh(2*h*kResults)^2 + 2509*cosh(2*h*kResults) - 1108)/(192*(cosh(2*h*kResults) - 1)^5);
    B53 = (57*cosh(2*h*kResults)^7 + 204*cosh(2*h*kResults)^6 - 53*cosh(2*h*kResults)^5 ...
        -782*cosh(2*h*kResults)^4 - 741*cosh(2*h*kResults)^3 - 52*cosh(2*h*kResults)^2 ...
        + 371*cosh(2*h*kResults) + 186)*9/((3*cosh(2*h*kResults) + 2)*(cosh(2*h*kResults) - 1)^6*128);

    switch modeNo
        case 1
            Results.a = a;
            Results.aw = a*(1+kResults^2.*a^2*B31+kResults^4.*a^4*B51);
            Hout = 2*a+2*(B31+B33)*kResults^2*a^3+2*(B51+B53+B55)*kResults^4*a^5;
            Results.H = Hout;
            omega2 = kResults^2*a^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
            omega4 = (a^4*kResults^4*(20*alpha1^5 + 112*alpha1^4 ...
        - 100*alpha1^3 - 68*alpha1^2 - 211*alpha1 + 328))/(32*(alpha1 - 1)^5);
            omega = omega0*(1+omega2+omega4);
            Results.omega = double(omega);
            Results.T = 2*pi/omega;
        case 2
            Results.aw = a; % in this case, users gave aw actually
            syms a0
            a0Fun = a - a0*(1+k^2.*a0^2*B31+k^4.*a0^4*B51);
            Results.a = double(vpasolve(a0Fun, a0, [0 inf]));
            a = Results.a;
            omega2 = kResults^2*a^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
            omega4 = (a^4*kResults^4*(20*alpha1^5 + 112*alpha1^4 ...
                - 100*alpha1^3 - 68*alpha1^2 - 211*alpha1 + 328))/(32*(alpha1 - 1)^5);
            omega = omega0*(1+omega2+omega4);
            Results.omega = double(omega);
            Results.T = 2*pi/omega;
            
            Hout = 2*a+2*(B31+B33)*kResults^2*a^3+2*(B51+B53+B55)*kResults^4*a^5;
            Results.H = Hout;
            
    end
elseif and(exist('H', 'var') ==1, exist('omega0', 'var') ==1)
    Results.omega0 = omega0;
    kResults = double(vpasolve(omega0^2 - g*k*tanh(k*h), k, [0 inf]));
    alpha1 = cosh(2*kResults*h);
    syms a 
    sigma = tanh(kResults*h);
    B33 = (27-9*sigma^2+9*sigma^4-3*sigma^6)/64/sigma^6;     
    B55 = (300*cosh(2*h*kResults)^8 + 1579*cosh(2*h*kResults)^7 + 3176*cosh(2*h*kResults)^6 + 2949*cosh(2*h*kResults)^5 ...
        + 1188*cosh(2*h*kResults)^4 + 675*cosh(2*h*kResults)^3 + 1326*cosh(2*h*kResults)^2 + 827*cosh(2*h*kResults) ...
        + 130)*5/((cosh(2*h*kResults) - 1)^6*(12*cosh(2*h*kResults)^2 + 11*cosh(2*h*kResults) + 2)*384);
    B31 = (3+8*sigma^2-9*sigma^4)/16/sigma^4;
    B51 = (121*cosh(2*h*kResults)^5 + 263*cosh(2*h*kResults)^4 +376*cosh(2*h*kResults)^3 ...
        - 1999*cosh(2*h*kResults)^2 + 2509*cosh(2*h*kResults) - 1108)/(192*(cosh(2*h*kResults) - 1)^5);
    B53 = (57*cosh(2*h*kResults)^7 + 204*cosh(2*h*kResults)^6 - 53*cosh(2*h*kResults)^5 ...
        -782*cosh(2*h*kResults)^4 - 741*cosh(2*h*kResults)^3 - 52*cosh(2*h*kResults)^2 ...
        + 371*cosh(2*h*kResults) + 186)*9/((3*cosh(2*h*kResults) + 2)*(cosh(2*h*kResults) - 1)^6*128);              
    Results.H = H;
    Hfun1 = H -(2*a+2*(B31+B33)*kResults^2*a^3+2*(B51+B53+B55)*kResults^4*a^5) ==0;
    aResults = double(vpasolve(Hfun1, [0 H]));
    omega2 = kResults^2*aResults^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
    omega4 = (aResults^4*kResults^4*(20*alpha1^5 + 112*alpha1^4 ...
        - 100*alpha1^3 - 68*alpha1^2 - 211*alpha1 + 328))/(32*(alpha1 - 1)^5);
    omega = omega0*(1+omega2+omega4);
    Results.a = aResults;
    Results.omega = omega;
    Results.T = 2*pi/omega;
    sigma = tanh(kResults*h);
        B31 = (3+8*sigma^2-9*sigma^4)/16/sigma^4;
        B51 = (121*cosh(2*h*kResults)^5 + 263*cosh(2*h*kResults)^4 +376*cosh(2*h*kResults)^3 ...
        - 1999*cosh(2*h*kResults)^2 + 2509*cosh(2*h*kResults) - 1108)/(192*(cosh(2*h*kResults) - 1)^5);
        Results.aw = aResults*(1+kResults^2.*aResults^2*B31+kResults^4.*aResults^4*B51);
  
else    
end

Results.k = kResults;
Results.L = 2*pi/kResults;
Results.c = Results.L/Results.T;
end
