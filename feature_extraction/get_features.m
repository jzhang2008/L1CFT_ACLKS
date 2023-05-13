function [feature_pixels, support_sz] = get_features(image, features, gparams, fg_size )

if ~ iscell(features)
    features = {features};
end;

[im_height, im_width, num_im_chan, num_images] = size(image);

colorImage = num_im_chan == 3;
num_features = length(features);

%compute total dimension of all features
tot_feature_dim = 0;
for n = 1:num_features
    
    if ~isfield(features{n}.fparams,'useForColor')
        features{n}.fparams.useForColor = true;
    end;
    
    if ~isfield(features{n}.fparams,'useForGray')
        features{n}.fparams.useForGray = true;
    end;
    
    if (features{n}.fparams.useForColor && colorImage) || (features{n}.fparams.useForGray && ~colorImage)
        tot_feature_dim = tot_feature_dim + features{n}.fparams.nDim;
    end;
    
end;

if nargin < 4 || isempty(fg_size)
    if gparams.cell_size == -1
        fg_size = size(features{1}.getFeature(image,features{1}.fparams,gparams));
    else
        fg_size = [floor(im_height/gparams.cell_size), floor(im_width/gparams.cell_size)];
    end
end

% temporary hack for fixing deep features
if gparams.cell_size == -1
    cf = features{1};
    if (cf.fparams.useForColor && colorImage) || (cf.fparams.useForGray && ~colorImage)
        [feature_pixels, support_sz] = cf.getFeature(image,cf.fparams,gparams);
    end;
else
    %compute the feature set
%     feature_pixels =cell(1,1,num_features);
%     ind =1;
%     currDim = 1;
%     for n=1:num_features
%         cf=features{n};
%         if (cf.fparams.useForColor && colorImage) || (cf.fparams.useForGray && ~colorImage)
%             num_blocks = 1;
%             feature_pixels{ind} = cf.getFeature(image,cf.fparams,gparams);
%             currDim = currDim + cf.fparams.nDim;
%         end;
%         ind = ind + num_blocks;
%     end
%     % Do feature normalization
%     if ~isempty(gparams.normalize_power) && gparams.normalize_power > 0
%         feature_pixels = cellfun(@(x) bsxfun(@times, x, ...
%         ((size(x,1)*size(x,2))^gparams.normalize_size * size(x,3)^gparams.normalize_dim ./ ...
%         (sum(abs(reshape(x, [], 1, 1, size(x,4))).^gparams.normalize_power, 1) + eps)).^(1/gparams.normalize_power)), ...
%             feature_pixels, 'uniformoutput', false);
%     end
%     feature_pixels=cell2mat(feature_pixels);

    feature_pixels = zeros(fg_size(1),fg_size(2),tot_feature_dim, num_images, 'single');
    currDim = 1;
    for n = 1:num_features
        cf = features{n};
        if (cf.fparams.useForColor && colorImage) || (cf.fparams.useForGray && ~colorImage)
            feature_pixels(:,:,currDim:(currDim+cf.fparams.nDim-1),:) = cf.getFeature(image,cf.fparams,gparams);
            currDim = currDim + cf.fparams.nDim;
        end;
    end;
    support_sz = [im_height, im_width];
end

end