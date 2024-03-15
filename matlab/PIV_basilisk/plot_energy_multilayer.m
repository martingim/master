function [] = plot_energy_multilayer(nx, nl)
%plots the total energy from the basilisk multilayer results

filename = sprintf('basilisk_results/energy_nx%d_nl%d.csv', nx, nl);
t = readtable(filename);
array  = table2array(t);
array(:,1) = array(:,1);
array(:,2) = array(:,2);
t = array(:,3);
% plot(t,array(:,1), 'DisplayName', 'ke');
% plot(t,array(:,2), 'DisplayName', 'gpe');
plot(t,(array(:,1) + array(:,2))/abs(array(1,1)+array(1,2)), 'DisplayName', sprintf('nx:%d, nl:%d', nx, nl));
legend();
xlabel('t');
ylabel('e/abs(e_0)'); 
title('ke + gpe');
hold on
end