function [crest_idx,v_mean, v_mean_lowpass] = find_crest(V)
%find_crest finds the crest of the wave by finding the x value with the
%lowest mean vertical velocity

V(isnan(V))=0;
v_abs = abs(V);
v_mean = mean(v_abs, 1);

%perform lowpass to smooth the data 0.001
v_mean_lowpass = lowpass(v_mean, 0.2);
%remove the leftmost and rightmost part of the data
v_mean_lowpass(1:floor(size(v_mean_lowpass,2)*0.25)) = NaN;
v_mean_lowpass(floor(size(v_mean_lowpass,2)*0.75):end) = NaN;

[ ~ , crest_idx] = min(v_mean_lowpass);

end