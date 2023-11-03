function [y_scaled, u_crest_scaled, alpha, crest_mask] = plot_velocities(run_number, pair_number, max_arrows)
%%plots the velocity fields

scale = 3;
height = 0.33;

%% load data and parameters
load('velocities.mat')
load('params.mat')

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
v = @(x, y) a*k*g/omega*exp(k*y).*sin(k*x-omega*t);


%% make New mask for quiver plots based on max_arows
quiver_idx = idx;
if size(quiver_idx, 2)>max_arrows || size(quiver_idx, 1)>max_arrows
    quiver_idx = zeros(size(quiver_idx));
    if size(quiver_idx, 1)>max_arrows
        interval_x = floor(size(quiver_idx, 1)/max_arrows);
    else
        interval_x = 1;
    end
    if size(quiver_idx, 2)>max_arrows
        interval_y = floor(size(quiver_idx, 2)/max_arrows);
    else
        interval_y = 1;
    end
    
    quiver_idx(1:interval_x:end, 1:interval_y:end) = 1;
    quiver_idx = logical(quiver_idx.*idx);
end


%% Quiver plots

%World
figure;
hold on;
quiver(xw(quiver_idx),yw(quiver_idx),Uw(quiver_idx),Vw(quiver_idx), scale);
legend('velocity')
title('world')
xlabel('x[m]')
ylabel('y[m]')
hold off;

%Analytical
figure;
quiver(xw(quiver_idx), yw(quiver_idx), u(xw(quiver_idx), yw(quiver_idx)), v(xw(quiver_idx), yw(quiver_idx)), scale)
title('analyctical')
xlabel('x[m]')
ylabel('y[m]')


%% Velocity profile

%find the crest
figure;
V = Vw;
V(isnan(Uw))=0;
v_abs = abs(V);
v_mean = mean(v_abs, 1);

NaN_sides = ones(size(v_mean));%vector to make the values to the left and right sides of the image NaN to find the internal minimum
for i = 1:size(v_mean, 2)
    if i<size(v_mean, 2)*0.25 || i>size(v_mean, 2)*0.75
        NaN_sides(i) = NaN;
    end
end

[ ~, crest_idx] = min(v_mean.*NaN_sides);
plot(xw(1,:), v_mean, 'r');
title('mean v over y')

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


%% Alpha plot

%calculate alpha
[M, N] = size(Uw);
V_norm = (Uw.^2 + Vw.^2).^0.5;
V_norm_mean = zeros(M,1)*NaN;

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
plot(alpha, yw(:,crest_idx)/height, 'x')
plot(exp(k*y_to_crest), y_to_crest/height)
legend('measured', 'theoretical')
title('Alpha plot')
xlabel('$\alpha$', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\frac{y}{h}$', 'interpreter', 'latex', 'FontSize', 20, 'rotation', 0)


%% Calculate sigma

sigma = 0;
m = 0;
for i=1:M
    m = m + max(idx(i,:)); % add one if there is a 1 in the mask on this row
    n = 1e-16;
    s = 0;
    for j=1:N
        if idx(i,j)==1
            n = n + 1;
            s = s + (V_norm(i,j)-V_norm_mean(i))^2;
        end
        
    end
    s = s/n;
    sigma = sigma + sqrt(s);
end
sigma = sigma/m

end
