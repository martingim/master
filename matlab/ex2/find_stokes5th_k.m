function [k] = find_stokes5th_k(run_number)
%%
image_names; %run this to get the frequencies
a = surface_height(run_number); %the height between the crest and the troughs
f = frequency(run_number);
H = 2*a;
T = 1/f;
Result = StokesDispSolver('h', 0.335, 'H', H,  'T', T, 'mode', 1);


try 
    load 'params.mat' params
catch
    params = containers.Map('KeyType', 'double', 'ValueType', 'any');
end

p = params(run_number);
p('k5') = Result.k;
p('omega5') = Result.omega;
params(run_number) = p;


save("params.mat", "params")

k = Result.k;

