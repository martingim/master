function [] = plot_alpha(run_number, pair_number)

load velocities.mat velocities
load params.mat params

%load parameters for analytical solution
p = params(run_number);
a = p('a');
k = p('k');
omega = p('omega');
%h = p('h');
h = 0.33;
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
[crest_idx, ~, ~] = find_crest(Vw);
%alpha analytical
alpha_analytical = @(y) exp(k*y);

%calculate alpha
[M, N] = size(Uw);
V_norm = (Uw.^2 + Vw.^2).^0.5;
V_norm_mean = zeros(M,1)*NaN;

%calculate the mean of V_norm for every y value
for i=1:M
    n = 0;
    s = 0;
    for j=1:N
        if idx(i,j)==1
            n = n + 1;
            s = s + V_norm(i,j);
        end
    end 
    V_norm_mean(i) = s/n;
end

alpha = omega/(a*g*k)*V_norm_mean;

figure
hold on
y_to_crest = yw(yw(:,1)<0.02, crest_idx);
plot(alpha, yw(:,crest_idx)/h, 'x')
plot(alpha_analytical(y_to_crest), y_to_crest/h)
legend('measured', 'theoretical')
title('Alpha plot')
xlabel('$\alpha$', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\frac{y}{h}$', 'interpreter', 'latex', 'FontSize', 20, 'rotation', 0)


