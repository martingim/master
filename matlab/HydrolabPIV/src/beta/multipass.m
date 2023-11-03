function [U,V,piv] = multipass(im1,mask1,im2,mask2)

opt = setpivopt('range',[-8 8 -8 8]);
piv = normalpass([],im1,mask1,im2,mask2,opt);

opt = setpivopt('range',[-3 3 -3 3]);  

piv = normalpass(piv,im1,mask1,im2,mask2,opt);
for i=1:2
    piv = distortedpass(piv,im1,mask1,im2,mask2,opt);
end

[U,V] = replaceoutliers(piv,0);