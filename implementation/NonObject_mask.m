function [Nonobj_m,im] = NonObject_mask(im,pos,target_sz)
ratio =0.75;
mask  = ones(size(im));
hs    = round(pos(1)-0.5*ratio*target_sz(1)):round(pos(1)+0.5*ratio*target_sz(1));
ws    = round(pos(2)-0.5*ratio*target_sz(2)):round(pos(2)+0.5*ratio*target_sz(2));
%check for out-of-bounds coordinates, and set them to the values at
%the borders
hs(hs<1)=0;
hs(hs>size(im,1))=0;
ws(ws<1)=0;
ws(ws>size(im,2))=0;
hs(hs==0)= [];
ws(ws==0)= [];
mask(hs,ws,:)= 0;
im(hs,ws,:)  = 0;
Nonobj_m = mask;
end