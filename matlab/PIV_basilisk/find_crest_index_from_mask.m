function [crest_idx] = find_crest_index_from_mask(y,mask)
%finds the index of the crest based on the mask and the world 
%coordinates

% disregards the leftmost and rightmost 20%
x_size = size(y, 2);
min_x_i = round(x_size*0.2);
max_x_i = round(x_size*0.8);
y(:,1:min_x_i)=-1;
y(:,max_x_i:end)=-1;

%finding the crest   +1 to avoid 0 being the max
[~, crest_idx] = max(max((y+1).*mask));


end