function [colorname_features]=get_colorname( im, fparam, gparam)
% extract color name features 
[im_height, im_width, num_im_chan, num_images] = size(im);
single_im = single(im);
% single_im=imresize(im,[im_height/gparam.cell_size, im_width/gparam.cell_size]);
% t_colorname_features=zeros(im_height/gparam.cell_size, im_width/gparam.cell_size, size(fparam.w2c,2), num_images,'single');
t_colorname_features=zeros(im_height, im_width, size(fparam.w2c,2), num_images,'single');
for k=1:num_images
    t_colorname_features(:,:,:,k)=im2c(single(single_im), fparam.w2c, -2);
end
if gparam.cell_size > 1
        colorname_features = average_feature_region(t_colorname_features,gparam.cell_size);
else
        colorname_features = t_colorname_features;
end;