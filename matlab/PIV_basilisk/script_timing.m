clear all
expression = '\d+\.?\d* CPU, [+-]?\d+\.?\d*( real,) [+-]?\d+\.?\d*e?[+-]?\d*';
%%
logfiles = [];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL7_r1_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL8_r1_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL9_r1_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL10_r1_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r1_log.txt"];
%% plot pointsteps/s for different LEVELS
LEVELS = [7 8 9 10 11 12];
logfiles = [];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL7_r11_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL8_r11_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL9_r11_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL10_r11_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL12_r11_log.txt"];



runtime = [];
cputime = [];
pointsteps = [];
% text = "# Quadtree, 4134 steps, 1242.42 CPU, 627.3 real, 2.66e+05 points.step/s, 42 var";


for logfile = logfiles
    text = readlines(logfile);
    [tokens, matches] = regexp(text,expression,'tokens',"match");
    for i=1:size(matches, 1)

        match = matches{i};
        if size(match,1)>0
            match
            times = sscanf(match, '%f CPU, %f real, %f');
            runtime = [runtime times(2)];
            cputime = [cputime times(1)];
            pointsteps = [pointsteps times(3)];
        end
    end
end

ylabel('seconds')
plot(LEVELS, pointsteps, '-s', 'DisplayName',"runtime 25s simulation")
xlabel('refinement LEVEL')
ylabel('points*steps/second')
ylim([0 max(pointsteps)])
legend()
hold on

runtime
cputime
pointsteps
%%
n = [1 2 3 4 5 6 7 8];
logfiles = [];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n1_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n2_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n3_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n4_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n5_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n6_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n7_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n8_log.txt"];


runtime = [];
cputime = [];
pointsteps = [];
% text = "# Quadtree, 4134 steps, 1242.42 CPU, 627.3 real, 2.66e+05 points.step/s, 42 var";


for logfile = logfiles
    text = readlines(logfile);
    [tokens, matches] = regexp(text,expression,'tokens',"match");
    for i=1:size(matches, 1)

        match = matches{i};
        if size(match,1)>0
            match
            times = sscanf(match, '%f CPU, %f real, %f');
            runtime = [runtime times(2)];
            cputime = [cputime times(1)];
            pointsteps = [pointsteps times(3)];
        end
    end
end


plot(n, pointsteps, '-s', 'DisplayName',"runtime 25s simulation")
xlabel('number of threads')
ylabel('points*steps/second')
ylim([0 max(pointsteps)])
title("Number of steps*points per second on machine with 8 cores")
legend()
hold on

runtime
cputime
pointsteps

%% PLOT NS from abacus
close all
n_threads = [1:14];
figure
hold on
for L = 10:12
    runtime = [];
    cputime = [];
    pointsteps = [];
    for n = n_threads
        text = readlines(sprintf("/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/results/LEVEL%d_r1_n%dlog.txt",L, n));                                                                                                     
        [tokens, matches] = regexp(text,expression,'tokens',"match");
        for i=1:size(matches, 1)
    
            match = matches{i};
            if size(match,1)>0
                match
                times = sscanf(match, '%f CPU, %f real, %f');
                runtime = [runtime times(2)];
                cputime = [cputime times(1)];
                pointsteps = [pointsteps times(3)];
            end
        end
    end
    plot(n_threads, runtime,'--', 'DisplayName',sprintf("Intel® Xeon® CPU E5-2695 NS LEVEL %d", L))
end


xlabel('number of threads')
ylabel('runtime [s]')
title('runtime for 5s simulation')
legend()


%% multilayer runtimes from abacus
% figure
hold on
n_threads = 1:14;
for L=10:12
    runtime = [];
    cputime = [];
    pointsteps = [];

    for n = n_threads
        text = readlines(sprintf("/home/martin/Documents/master/basilisk/2d_piston/multilayer-piston/results/LEVEL%d_nl_20_r1_n%d_log.txt", L, n));
                                                                                                             
        [tokens, matches] = regexp(text,expression,'tokens',"match");
        for i=1:size(matches, 1)
    
            match = matches{i};
            if size(match,1)>0
                match
                times = sscanf(match, '%f CPU, %f real, %f');
                runtime = [runtime times(2)];
                cputime = [cputime times(1)];
                pointsteps = [pointsteps times(3)];
            end
        end
    end
    
    
    
    plot(n_threads, runtime,'DisplayName',sprintf("Intel® Xeon® CPU E5-2695 multilayer LEVEL %d", L))
end
xlabel('number of threads')
ylabel('runtime [s]')
title('runtime for 5s simulation on cpu with 18 cores')
legend()
fontsize(20, "points")
% runtime
% cputime
% pointsteps


%% multilayer runtimes my pc
figure
hold on
n_threads = 1:8;
for L=10:12
    runtime = [];
    cputime = [];
    pointsteps = [];

    for n = n_threads
        text = readlines(sprintf("/home/martin/Documents/master/basilisk/2d_piston/multilayer-piston/LEVEL%d_nl_20_r1_n%d_log.txt", L, n));
                                                                                                             
        [tokens, matches] = regexp(text,expression,'tokens',"match");
        for i=1:size(matches, 1)
    
            match = matches{i};
            if size(match,1)>0
                match
                times = sscanf(match, '%f CPU, %f real, %f');
                runtime = [runtime times(2)];
                cputime = [cputime times(1)];
                pointsteps = [pointsteps times(3)];
            end
        end
    end
    
    
    
    plot(n_threads, runtime,'DisplayName',sprintf("Intel® Core™ i7-9700K multilayer LEVEL %d", L))
end
xlabel('number of threads')
ylabel('runtime [s]')
title('runtime for 5s simulation on machine with 8 cores')
legend()
fontsize(20, "points")
% runtime
% cputime
% pointsteps

%% PLOT NS from my pc

n_threads = [1:8];

for L = 10:12
    runtime = [];
    cputime = [];
    pointsteps = [];
    for n = n_threads
        text = readlines(sprintf("/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL%d_r1_n%d_log_5s.txt",L, n));
        [tokens, matches] = regexp(text,expression,'tokens',"match");
        for i=1:size(matches, 1)
    
            match = matches{i};
            if size(match,1)>0
                match
                times = sscanf(match, '%f CPU, %f real, %f');
                runtime = [runtime times(2)];
                cputime = [cputime times(1)];
                pointsteps = [pointsteps times(3)];
            end
        end
    end
    plot(n_threads, runtime,'DisplayName',sprintf("Intel® Core™ i7-9700K NS LEVEL %d", L))
end


xlabel('number of threads')
ylabel('runtime [s]')
title('runtime for 5s simulation')
legend()
%%
fontsize(20, "points")