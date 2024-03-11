function [mask] = mask_image(im1)
%makes a transform from pixel to coordinate
%   Detailed explanation goes here
% create a polygonal mask

figure;
set(gca,'Ydir','normal')
imagesc(im1);
h = impoly();
mask = h.createMask();
end