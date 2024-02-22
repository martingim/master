function [k] = find_stokes5th_k(run_number, image_params)
%%


try 
    load 'params.mat' params
    p = params(run_number);
    a = p('a');
catch
    params = containers.Map('KeyType', 'double', 'ValueType', 'any');
    a = surface_height(run_number); %the height between the crest and the troughs
end

frequency = image_params('frequency');
water_depth = image_params('water_depth');
h = water_depth(run_number);
f = frequency(run_number);
H = 2*a;
T = 1/f;
Result = StokesDispSolver('h', h, 'H', H,  'T', T, 'mode', 1);





p('k5') = Result.k;
p('omega5') = Result.omega;
params(run_number) = p;


save("params.mat", "params")

k = Result.k;

