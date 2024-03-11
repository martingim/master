target_k = 7.14;
k = 0;
f = 1;
while abs(target_k - k)>0.0001
    f = f + 0.1*(target_k - k);
    T = 1/f;
    h = 0.6;
    a = 0.05;
    omega = 2*pi/T;
    g = 9.81;
    t = 0;
    H = a*2;
    Result = StokesDispSolver('h', h, 'H', H,  'T', T, 'mode', 1);
    
    k = Result.k;
end


z = 0:-0.001:-h;

[phi, u,~,~,~,~,~, w,~,~,~,~,~] = StokesU(Result.k, Result.h, Result.a, 0, z);


plot(u, z)
hold on 
[omega,T,k,lambda,Cp,Cg] = wparam(f, h);
u_old = @(x, y) a*k*g/omega*exp(k*y).*cos(k*x-omega*t);
v_old = @(x, y) a*k*g/omega*exp(k*y).*sin(k*x-omega*t);
plot(u_old(0,z), z)


legend('5th order', 'second')
