function [alpha, alpha_theoretical, y_at_crest, y_below_crest_scaled] = plot_alpha(run_number, pair_number, image_params, create_plot)

load velocities.mat velocities
load params.mat params

%load parameters for analytical solution
p = params(run_number);
a = p('a');
a_std = p('std_a');
k = p('k');
omega = p('omega');
water_depth = image_params('water_depth');
h = water_depth(run_number);
g = 9.81;
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
y_at_crest = yw(:,crest_idx)/h;
y_to_crest = yw(yw(:,1)<=a, crest_idx);
y_below_crest_scaled = y_to_crest/h;
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

%calculate the standard deviation
V_norm_std = zeros(M,1)*NaN;
for i=1:M
    n = 0;
    s = 0;
    for j=1:N
        if idx(i,j)==1
            n = n + 1;
            s = (V_norm_mean(i)-V_norm(i,j)).^2;
        end
    end 
    V_norm_std(i) = s.^.5/n;
end


alpha = omega/(a*g*k)*V_norm_mean;
alpha_std = omega/(a*g*k)*V_norm_std;
alpha_theoretical = alpha_analytical(y_to_crest);
if create_plot
    figure
    hold on
    plot(alpha, y_at_crest, 'x')
    plot(alpha_theoretical, y_below_crest_scaled)
    legend('measured', 'theoretical')
    title('Alpha plot')
    xlabel('$\alpha=\frac{\omega}{agk}(u^2+v^2)^{\frac{1}{2}}$', 'interpreter', 'latex', 'FontSize', 20)
    ylabel('$\frac{y}{h}$', 'interpreter', 'latex', 'FontSize', 20, 'rotation', 0)
end

end