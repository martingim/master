clear all
expression = '\d+\.?\d* CPU, [+-]?\d+\.?\d*( real,) [+-]?\d+\.?\d*e?[+-]?\d*';
LEVELS = [1 2 3 4 5 6 7 8];
%%
logfiles = [];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL7_r1_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL8_r1_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL9_r1_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL10_r1_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r1_log.txt"];
%%
logfiles = [];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL7_r11_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL8_r11_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL9_r11_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL10_r11_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL12_r11_log.txt"];
%%
logfiles = [];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n1_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n2_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n3_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n4_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n5_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n6_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n7_log.txt"];
logfiles = [logfiles "/home/martin/Documents/master/basilisk/2d_piston/boundary-piston/LEVEL11_r11_n8_log.txt"];

%%

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
plot(LEVELS, runtime, '-s', 'DisplayName',"runtime 25s simulation")
xlabel('maximum refinement level')
ylim([0 max(runtime)])
legend()
hold on

runtime
cputime
pointsteps