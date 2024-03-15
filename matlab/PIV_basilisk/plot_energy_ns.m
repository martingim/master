function [] = plot_energy_ns(N, n)
%plots the total energy from the basilisk multilayer results

filename = sprintf('basilisk_results/energy_ns_coarseN%d_fineN%d.csv', N, n);
t = readtable(filename);
array  = table2array(t);
array(:,1) = array(:,1);
array(:,2) = array(:,2);
t = array(:,3);
%plot(t,array(:,1), 'DisplayName', 'ke');
%plot(t,array(:,2), 'DisplayName', 'gpe');
plot(t,(array(:,1) + array(:,2))/abs(array(1,1)+array(1,2)), 'DisplayName', sprintf('N:%d, refinement:%d', N,n));
legend();
xlabel('t');
ylabel('e/abs(e_0)'); 
title('ke + gpe');
hold on
end