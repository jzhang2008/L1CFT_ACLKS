function results = L1CFT_ACLKS(params)
%% Initialization
project_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(project_dir, 'utils'));
addpath(fullfile(project_dir, 'implementation'));
addpath(fullfile(project_dir, 'feature_extraction'));
% Get sequence info
[seq, im] = get_sequence_info(params.seq);
params = rmfield(params, 'seq');
if isempty(im)
    seq.rect_position = [];
    [~, results] = get_sequence_results(seq);
    return;
end
% parameters
search_area_scale = params.search_area_scale;
% padding = params.padding;                         	%extra area surrounding the target
output_sigma_factor = params.output_sigma_factor;	%spatial bandwidth (proportional to target)
% lambda = params.lambda;
learning_rate = params.learning_rate;
nScales = params.number_of_scales;
scale_step = params.scale_step;
% scale_sigma_factor = params.scale_sigma_factor;
% scale_model_max_area = params.scale_model_max_area;
refinement_iterations = params.refinement_iterations;
filter_max_area = params.filter_max_area;
interpolate_response = params.interpolate_response;
% num_GS_iter = params.num_GS_iter;


% video_path = params.video_path;
img_files = params.img_files;
pos = floor(params.init_pos);                                              % Init potition
% context position
pos_left    = floor(params.init_pos_left);
pos_right    = floor(params.init_pos_right);
pos_bottom    = floor(params.init_pos_bottom);
pos_top    = floor(params.init_pos_top);

pos_lefttop =  params.init_pos_lefttop ;
pos_righttop = params.init_pos_righttop ;
pos_leftbottom = params.init_pos_leftbottom ;
pos_rightbottom = params.init_pos_rightbottom;

pos_context = [pos_left; pos_right; pos_bottom; pos_top; pos_lefttop; pos_righttop; pos_leftbottom; pos_rightbottom];
keystep = params.keystep;

target_sz = params.wsize;
params.init_sz = target_sz;

debug = params.debug;
% visualization = params.visualization || debug;
%Feature settings
features = params.t_features;
% Set default parameters
params = init_default_params(params);

% Global feature parameters
if isfield(params, 't_global')
    global_fparams = params.t_global;
else
    global_fparams = [];
end
global_fparams.use_gpu = params.use_gpu;
global_fparams.gpu_id = params.gpu_id;

% Define data types
if params.use_gpu
    params.data_type = zeros(1, 'single', 'gpuArray');
else
    params.data_type = zeros(1, 'single');
end
params.data_type_complex = complex(params.data_type);

global_fparams.data_type = params.data_type;
%visualization = params.visualization;
% Check if color image
im = imread(img_files{1});
if size(im,3) == 3
    if all(all(im(:,:,1) == im(:,:,2)))
        is_colorImage = false;
    else
        is_colorImage = true;
    end
else
    is_colorImage = false;
end
if size(im,3) > 1 && is_colorImage == false
    im = im(:,:,1);
end


init_target_sz = target_sz;
%set the feature ratio to the feature-cell size
featureRatio = global_fparams.cell_size;
search_area = prod(init_target_sz / featureRatio * search_area_scale);

% when the number of cells are small, choose a smaller cell size
if isfield(global_fparams, 'cell_selection_thresh')
    if search_area < global_fparams.cell_selection_thresh * filter_max_area
        params.t_global.cell_size = min(featureRatio, max(1, ceil(sqrt(prod(init_target_sz * search_area_scale)/(global_fparams.cell_selection_thresh * filter_max_area)))));
        
        featureRatio = global_fparams.cell_size;
        search_area = prod(init_target_sz / featureRatio * search_area_scale);
    end
end


% if search_area > filter_max_area
%     currentScaleFactor = sqrt(search_area / filter_max_area);
% else
%     currentScaleFactor = 1.0;
% end
% global_feat_params = params.t_global;
search_area = prod(init_target_sz * params.search_area_scale);
if search_area > params.max_image_sample_size
    currentScaleFactor = sqrt(search_area / params.max_image_sample_size);
elseif search_area < params.min_image_sample_size
    currentScaleFactor = sqrt(search_area / params.min_image_sample_size);
else
    currentScaleFactor = 1.0;
end
% target size at the initial scale
base_target_sz = target_sz/currentScaleFactor;
%window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        img_sample_sz = floor( base_target_sz * search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        img_sample_sz = repmat(sqrt(prod(base_target_sz * search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
    case 'fix_padding'
        img_sample_sz = base_target_sz + sqrt(prod(base_target_sz * search_area_scale) + (base_target_sz(1) - base_target_sz(2))/4) - sum(base_target_sz)/2; % const padding
    otherwise
        error('Unknown "params.search_area_shape". Must be ''proportional'', ''square'' or ''fix_padding''');
end
% set the size to exactly match the cell size
img_sample_sz = round(img_sample_sz / featureRatio) * featureRatio;
use_sz = floor(img_sample_sz/featureRatio);
% the real crop patch range
params.real_area_scale =img_sample_sz./base_target_sz;

[features, global_fparams, feature_info] = init_features(features, global_fparams, is_colorImage, img_sample_sz, 'exact');
% Get feature specific parameters
feature_extract_info = get_feature_extract_info(features);

% window size, taking padding into account
% sz = floor(base_target_sz * (1 + padding));

% desired translation filter output (gaussian shaped), bandwidth
% proportional to target size
% construct the label function
output_sigma = sqrt(prod(floor(base_target_sz/featureRatio))) * output_sigma_factor;
rg = circshift(-floor((use_sz(1)-1)/2):ceil((use_sz(1)-1)/2), [0 -floor((use_sz(1)-1)/2)]);
cg = circshift(-floor((use_sz(2)-1)/2):ceil((use_sz(2)-1)/2), [0 -floor((use_sz(2)-1)/2)]);
% [rs, cs] = ndgrid(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), -floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2));
[rs, cs] = ndgrid( rg,cg);
y = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf = single(fft2(y));

if interpolate_response == 1
    interp_sz = use_sz * featureRatio;
else
    interp_sz = use_sz;
end

% desired scale filter output (gaussian shaped), bandwidth proportional to
% number of scales
% scale_sigma = nScales/sqrt(33) * scale_sigma_factor;
% ss = (1:nScales) - ceil(nScales/2);
% ys = exp(-0.5 * (ss.^2) / scale_sigma^2);
% ysf = single(fft(ys));

% store pre-computed translation filter cosine window
cos_window = single(hann(use_sz(1)) * hann(use_sz(2))');
% the search area size
% support_sz = prod(use_sz);
% Calculate feature dimension

% store pre-computed scale filter cosine window
% if mod(nScales,2) == 0
%     scale_window = single(hann(nScales+1));
%     scale_window = scale_window(2:end);
% else
%     scale_window = single(hann(nScales));
% end;
% if params.use_reg_window
%     % create weight window
%     % normalization factor
%     reg_scale = 0.5 * base_target_sz/featureRatio;
%     % construct grid
%     wrg = -(use_sz(1)-1)/2:(use_sz(1)-1)/2;
%     wcg = -(use_sz(2)-1)/2:(use_sz(2)-1)/2;
%     [wrs, wcs] = ndgrid(wrg, wcg);
%     % construct the regukarization window
%     reg_window = exp(-(abs(wrs/reg_scale(1)).^params.reg_window_power + abs(wcs/reg_scale(2)).^params.reg_window_power)/12);
% end
% Define spatial regularization windows
if params.use_reg_window
    reg_scale = floor(1.3*base_target_sz/featureRatio);    
    reg_window = ones(use_sz) * 1e5;
    range = zeros(numel(reg_scale), 2);
    
    % determine the target center and range in the regularization windows
    for j = 1:numel(reg_scale)
        range(j,:) = [0, reg_scale(j) - 1] - floor(reg_scale(j) / 2);
    end
    center = floor((use_sz + 1)/ 2) + mod(use_sz + 1,2);
    range_h = (center(1)+ range(1,1)) : (center(1) + range(1,2));
    range_w = (center(2)+ range(2,1)) : (center(2) + range(2,2));
    
    reg_window(range_h, range_w) = 1e-10;
end
params.reg_window = reg_window;
small_filter_sz=floor(base_target_sz/featureRatio);
[sx,sy,~] = get_subwindow_no_window(reg_window, floor(use_sz/2) , small_filter_sz);
% reg_window(sx,sy)=ones(small_filter_sz);
% scale factors
if nScales > 0
    scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2));
    
    scaleFactors = scale_step .^ scale_exp;
    
    %force reasonable scale changes
    min_scale_factor = scale_step ^ ceil(log(max(5 ./ img_sample_sz)) / log(scale_step));
    max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
end

% compute the resize dimensions used for feature extraction in the scale
% estimation
if interpolate_response >= 3
    % Pre-computes the grid that is used for socre optimization
%     ky = -floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2);
%     kx = -floor((use_sz(2) - 1)/2): ceil((use_sz(2) - 1)/2);
%     kx = kx';
    ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
    kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';
    newton_iterations = params.newton_iterations;
end
% currentScaleFactor = 1;

% to calculate precision

% to calculate FPS
seq.time = 0;
% the index of dictionary
dict_idx = ones(1,2);
% Define the learning variables
xl_context = cell(1, 8);
xlw_context = cell(1, 8);
xlf_context = cell(1, 8);
S_xx_context = cell(1, 8);
Dict         = {cell(1,10),cell(1,10)};
model_xf     = cell(1, 2);
cf_f         = cell(1, 2);
g_pre_f_new      = cell(1, 2);
% Allocate
scores_fs_feat = cell(1,2);
buffer_st_xl   = cell(1,15);
buffer_gf      = cell(1,15);
buffer_idx     =0;
% find maximum and minimum scales
% im = imread([video_path img_files{1}]);
% min_scale_factor = scale_step ^ ceil(log(max(5 ./ sz)) / log(scale_step));
% max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
% multires_pixel_template = zeros(img_sample_sz(1), img_sample_sz(2), size(im,3), nScales, 'uint8');
while true
    %load image
%     ss=[params.video_path,img_files{frame}];
%     im = imread(img_files{frame});
%     if size(im,3) > 1 && colorImage == false
%         im = im(:,:,1);
%     end
    if seq.frame>0
        [seq, im] = get_sequence_frame(seq);
        if isempty(im)
            break;
        end
        if size(im,3) > 1 && is_colorImage == false
            im = im(:,:,1);
        end
    else
        seq.frame=1;
    end
    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Target localization step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Do not estimate translation and scaling on the first frame, since we 
    % just want to initialize the tracker there
    if seq.frame > 1
        old_pos = inf(size(pos));
        iter = 1;
        %translation search
        while iter <= refinement_iterations && any(old_pos ~= pos)
            sample_pos = round(pos);
            sample_scale = currentScaleFactor*scaleFactors;
            xt = extract_features(im, sample_pos, sample_scale, features, global_fparams, feature_extract_info);                     
            % Do windowing of features
            xtw = cellfun(@(feat_map) bsxfun(@times, feat_map, cos_window), xt, 'uniformoutput', false);
            
            % Compute the fourier series
            xtf = cellfun(@fft2, xtw, 'uniformoutput', false);
            responsef=cellfun(@(hf,x) permute(sum(bsxfun(@times,conj(hf), x), 3),[1 2 4 3]),cf_f,xtf,'uniformoutput', false);
            if interpolate_response == 2
                % use dynamic interp size
                interp_sz = floor(size(y) * featureRatio * currentScaleFactor);
            end
            %--------------------------------------------------------------
%             % not consider the adaptive fusion by selecting the below several lines codes
            responsef{1} = responsef{1}+responsef{2};
            responsef{2} =[];
            responsef(cellfun(@isempty,responsef))=[];
            %-----------------------------------------------------------------------------
            responsef_padded =cellfun(@(resp) resizeDFT2(resp,interp_sz),responsef,'uniformoutput', false);
            response         =cellfun(@(resp) ifft2(resp,'symmetric'),responsef_padded,'uniformoutput', false);
            [row,col,sind]   =cellfun(@(resp) ind2sub(size(resp), find(resp == max(resp(:)), 1)),response,'uniformoutput', false);
             if length(sind)>1
                for k=1:length(sind)
                     PSR(k)=CalPSR(response{1,k}(:,:,sind{k}),11,use_sz);
%                    PSR(k)=CalAPCE(response{1,k}(:,:,sind{k}));
                end
                
%                 if min(PSR)<10
%                     learning_rates=0;
%                 else
%                     learning_rates=params.learning_rate;
%                 end
                scores_fs_padded   =PSR(1)/sum(PSR)*responsef_padded{1}(:,:,sind{1})+PSR(2)/sum(PSR)*responsef_padded{2}(:,:,sind{2});
                scores_fs          =PSR(1)/sum(PSR)*response{1}(:,:,sind{1})+PSR(2)/sum(PSR)*response{2}(:,:,sind{2});
                scale_change_factor=PSR(1)/sum(PSR)*scaleFactors(sind{1})+PSR(2)/sum(PSR)*scaleFactors(sind{2});
            else
                scores_fs_padded   = responsef_padded{1}(:,:,sind{1});
                scores_fs          = response{1}(:,:,sind{1});
                scale_change_factor=scaleFactors(sind{1});
            end
            % find maximum
            [disp_row, disp_col, ~] = resp_newton(scores_fs, scores_fs_padded, newton_iterations, ky, kx, use_sz);
            
            % Compute the translation vector in pixel-coordinates and round
            % to the closest integer pixel.
            translation_vec = [disp_row, disp_col] * featureRatio * currentScaleFactor *  scale_change_factor;
            % set the scale
            currentScaleFactor = currentScaleFactor * scale_change_factor;
            % set the scale
%             currentScaleFactor = currentScaleFactor * scaleFactors(sind);
           
            % adjust to make sure we are not to large or to small
            if currentScaleFactor < min_scale_factor
                currentScaleFactor = min_scale_factor;
            elseif currentScaleFactor > max_scale_factor
                currentScaleFactor = max_scale_factor;
            end
            
            % update position
            old_pos = pos;
            pos = sample_pos + translation_vec;
            
            iter = iter + 1;
        end
    end
    if seq.frame>1
        inv_trans_vec=[-disp_row,-disp_col];
        [pos_local_trans,params.distractor_num,local_max]=generate_distractor_pos(scores_fs,small_filter_sz,inv_trans_vec,params);
        local_trans_vec = pos_local_trans* featureRatio * currentScaleFactor;
        local_pos       = ones(params.distractor_num,1)*pos+local_trans_vec;
        target_sz       = floor(base_target_sz * currentScaleFactor);
        dist_w          = 1-sqrt(sum(local_trans_vec.^2,2))/sqrt(sum(target_sz.^2));
        params.weight   = dist_w'.*local_max*params.weight_const; %
%         figure(2),imshow(im);
%         hold on;
%         rect_pos  = [pos([2,1]) - (target_sz([2,1]) - 1)/2, target_sz([2,1])];
%         rectangle('Position',rect_pos, 'EdgeColor','g', 'LineWidth',2);
%         for i=1:params.distractor_num
%             local_rect_pos(i,:)  = [local_pos(i,[2,1]) - (target_sz([2,1]) - 1)/2, target_sz([2,1])];
%             rectangle('Position',local_rect_pos(i,:), 'EdgeColor','b', 'LineWidth',2);
%         end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Model update step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract the training sample feature map for the translation filter
%     xl = get_translation_sample(im, pos, sz, currentScaleFactor, cos_window);
    sample_pos = round(pos);
    xl = extract_features(im, sample_pos, currentScaleFactor, features, global_fparams, feature_extract_info);
    % do windowing of features
    xlw = cellfun(@(feat_map) bsxfun(@times, feat_map, cos_window), xl,  'uniformoutput', false);
    % compute the fourier series
    xlf = cellfun(@fft2, xlw, 'uniformoutput', false);
if params.adaptive_select
%     % if exist the distractor, crop the corresponding patch
    if seq.frame==1
        params.distractor_num=0;
    end
    if (seq.frame>1)
        if params.distractor_num>0
           for p=1:params.distractor_num
            % extract features
              xl_context{p} = extract_features(im, local_pos(p,:), currentScaleFactor, features, global_fparams, feature_extract_info);
            % do windowing of features
              xlw_context{p} =  cellfun(@(feat_map) bsxfun(@times, feat_map, cos_window), xl_context{p}, 'uniformoutput', false);
            % compute the fourier series
              xlf_context{p} = cellfun(@fft2, xlw_context{p}, 'uniformoutput', false);
           end
        end
    end
else
%     %%%% keyframe context learning
    if ((mod(seq.frame,keystep * params.Period) == 0)||(seq.frame==1))
        for p=1:8
            % extract features
            xl_context{p} = extract_features(im, pos_context(p,:), currentScaleFactor, features, global_fparams, feature_extract_info);
            % do windowing of features
            xlw_context{p} =  cellfun(@(feat_map) bsxfun(@times, feat_map, cos_window), xl_context{p}, 'uniformoutput', false);
            % compute the fourier series
            xlf_context{p} = cellfun(@fft2, xlw_context{p}, 'uniformoutput', false);
        end
    end
end
%   % select the optimal 
    if seq.frame==2
       g_pre_f_new = g_pre_f;
    end
    if seq.frame>2
       buffer_xl      =xl{1}(range_h, range_w,:);
       vec_xl         =buffer_xl(:);
       buffer_st_xl(cellfun(@isempty,buffer_st_xl))=[];
       buffer_fullM=cell2mat(buffer_st_xl);
       corrValue   =corr(vec_xl,buffer_fullM);                             % calculate the corrlation coefficients between the current patch and each buffer vector
       [~,corr_idx]=max(corrValue);
%        g_pre_f_new =buffer_gf{corr_idx};
       if length(corrValue)>=15
           [~,min_corr_idx] = min(corrValue(2:end));
       end
    end
    if seq.frame>1
      g_pre_f_new = g_pre_f;
    end
    % train the CF model for each feature
    [cf_f,g_pre_f]=train_filter(yf,xlf,xlf_context,cf_f,g_pre_f_new,params,seq,use_sz);
    %keyfilter update
%     if(mod(seq.frame,params.keystep)==0||seq.frame==1)
%        g_pre_f_new = g_pre_f;
%     end
    %adaptive select the dynamic temporal regulirization
    if buffer_idx<15
        buffer_idx  = buffer_idx +1;
        buffer_xl  = xl{1}(range_h, range_w,:);
        buffer_st_xl{buffer_idx}=buffer_xl(:);
        buffer_gf{buffer_idx}   =g_pre_f;
    else
%         buffer_M   =cell2mat(buffer_st_xl(2:15));
%         Gram_xl    =corr(buffer_M);
%         buffer_xt_new = buffer_st_xl;
%         buffer_xl  = xl{1}(range_h, range_w,:);
%         buffer_xt_new{2}=buffer_xl(:);
%         buffer_M_new  = cell2mat(buffer_xt_new(2:15));
%         Gram_xl_new    =corr(buffer_M_new);
%         dGram_xl       = det(Gram_xl);
%         dGram_xl_new   = det(Gram_xl_new);
%         if dGram_xl_new>dGram_xl
%            buffer_st_xl{2} = [];
%            buffer_st_xl(cellfun(@isempty,buffer_st_xl))=[];
%            buffer_xl  = xl{1}(range_h, range_w,:);
%            buffer_st_xl{15}= buffer_xl(:);
%            buffer_gf{2} = [];
%            buffer_gf(cellfun(@isempty,buffer_gf))=[];
%            buffer_gf{15}= g_pre_f;
%         end
        buffer_gf{min_corr_idx+1} = [];
        buffer_gf(cellfun(@isempty,buffer_gf))=[];
        buffer_gf{15}= g_pre_f;

    end
    
%   for k= 1: numel(xlf)
%       if(seq.frame==1)
%           model_xf{k} = xlf{k};
%       else
%           model_xf{k} =(1-learning_rates)*model_xf{k}+learning_rates*xlf{k};
%       end
%       [cf_f{k},g_pre_f{k}]=train_filter(yf,xlf{k},xlf_context,cf_f{k},g_pre_f{k},params,seq);
%       model_xf{k} = xlf{k};
%       % filter update of every Layer model
%       g_f        = single(zeros(size(xlf{k})));
%       h          = zeros(size(xlf{k}));
%       h_f        = g_f;
%       l_f        = g_f;
%       v          = zeros(1,length(Dict{k}))';  %Auxiliary Variables
%       sx         = zeros(1,length(Dict{k}))';  %real variables
%       eta        = zeros(1,length(Dict{k}))';
%       mu         = 1;
%       alpha      = 1.1;
%       mumax      = 10;
%       omega      = reg_window;
%       gamma      = 10;
%       theta      = 1;
%       lambda     = 0.05;
%       RELTOL     = 0.5*1*1e-2;
%       Iter       = 50;%100;
%       i          = 1;
%       res        = zeros(1,Iter);
%       S_xx0 = sum(conj(model_xf{k}) .* model_xf{k}, 3);
%       S_xy = bsxfun(@times, yf, conj(model_xf{k}));
%     %keyframe context learning
%     if ((mod(seq.frame,keystep * params.Period) == 0)||(seq.frame==1))
%         for p=1:8
%             S_xx_context{p}=sum(conj(xlf_context{p}{1,k}).* xlf_context{p}{1,k},3);
%         end
%         S_xx = S_xx0 + params.yta1*(S_xx_context{1} + S_xx_context{2})...
%              + params.yta2 * (S_xx_context{3} + S_xx_context{4})...
%              + params.yta3 * (S_xx_context{5} + S_xx_context{6} + S_xx_context{7} + S_xx_context{8});
%     else
%        S_xx = S_xx0;  
%     end
%     %ADMM
%     while(i<Iter)
%         %   solve for G- please refer to the paper for more details
%            g_f = bsxfun(@times,(S_xy-l_f+mu*h_f),1./(S_xx+mu));
%         %   solve for H
% %            h   = shrinkage(real(ifft2(g_f)+ifft2(l_f)/mu),omega/mu);
%         if seq.frame<=10
%            h   = shrinkage(real(ifft2(g_f)+ifft2(l_f)/mu),omega/mu);
%         else
%            Q   = cell_DotProduct_vec(Dict{k},v'); 
%            h   = shrinkage(real((mu*ifft2(g_f)+ifft2(l_f)+gamma*Q)/(mu+gamma)),omega/(mu+gamma));
%            [hnum,hden] = calTensorProduct(Dict{k},h);
% %            v   = (gamma*hden+theta*eye(10))\(gamma*hnum-eta+theta*sx);
%            v   = (gamma*hden+theta*eye(10))\(gamma*hnum+eta+theta*sx);
% %           sx  = v>=max(v);
% %            sx  = shrinkage(v-eta/theta,lambda/theta);
%             sx  = max(v-(eta+lambda)/theta,0);
%            eta = eta+theta*(sx-v);
%            theta = min(1.5*theta,100);
%         end
%          h_f = fft2(h);
%          % update L
%          l_f =l_f+mu*(g_f-h_f);
%          if max(max(sum( abs(ifft2(g_f)-h),2)))<RELTOL||i>=Iter  %&& max(max(sum( abs(ifft2(f_w-f_w_old)),2)))<RELTOL)
%              break;
%          end
%          res(i)=max(max(sum( abs(ifft2(g_f)-h),2)));
%         %   update mu- betha = 10.
% %         mu = min(alpha * mu, mumax);
%          i=i+1;
%     end
%   
%     if seq.frame<=10
%         Dict{k}{dict_idx(k)}= h;
%         dict_idx(k)      = dict_idx(k)+1;
%     else
%         TDict          =Dict{k};
%         TDict{1}       = [];
%         TDict(cellfun(@isempty,TDict))=[];
%         
%         TDict{10}      = h; 
%         Dict{k}        =TDict;
%     end
%     cf_f{k}  = h_f;
%   end
    cf_f(cellfun(@isempty,cf_f))=[];
    % calculate the new target size
    target_sz = floor(base_target_sz * currentScaleFactor);
    % update the position of context patches
    pos_left = pos - [target_sz(1,2)/2,0];
    pos_right = pos + [target_sz(1,2)/2,0];
    pos_bottom = pos - [0, target_sz(1,1)/2];
    pos_top = pos + [0, target_sz(1,1)/2];

    pos_lefttop = pos - [target_sz(1,2)/2,-1*target_sz(1,1)/2];
    pos_righttop = pos + [target_sz(1,2)/2,target_sz(1,1)/2];
    pos_leftbottom = pos - [target_sz(1,2)/2, target_sz(1,1)/2];
    pos_rightbottom = pos + [-1*target_sz(1,2)/2,target_sz(1,1)/2];
    
    pos_context = [pos_left; pos_right; pos_bottom; pos_top; pos_lefttop; pos_righttop; pos_leftbottom; pos_rightbottom];
    %save position and calculate FPS
    tracking_result.center_pos = double(pos);
    tracking_result.target_size = double(target_sz);
    seq = report_tracking_result(seq, tracking_result);
    seq.time = seq.time + toc();
    %save position
%     positions(frame,:) = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
%     
%     time = time + toc;
    
    
    %visualization
%     if visualization ==1
%         rect_position_vis = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
%         if isToolboxAvailable('Computer Vision Toolbox')
%             im = insertShape(im, 'Rectangle', rect_position_vis, 'LineWidth', 4, 'Color', 'red');
% %                 im = insertShape(im, 'Rectangle', rect_position_padded, 'LineWidth', 4, 'Color', 'yellow');
%                 % Display the annotated video frame using the video player object. 
%             step(params.videoPlayer, im);
%         end
%     end
     %% Visualization
    if params.visualization
        rect_position_vis = [pos([2,1]) - (target_sz([2,1]) - 1)/2, target_sz([2,1])];
        outrect_position  = [pos([2,1]) - params.real_area_scale([2,1]).*(target_sz([2,1]) - 1)/2, params.real_area_scale([2,1]).*target_sz([2,1])];
        im_to_show = double(im)/255;
        if size(im_to_show,3) == 1
            im_to_show = repmat(im_to_show, [1 1 3]);
        end      
        imagesc(im_to_show);
        hold on;
        rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
        rectangle('Position',outrect_position, 'EdgeColor','r', 'LineWidth',3);
        text(10, 10, [int2str(seq.frame) '/'  int2str(numel(img_files))], 'color', [0 1 1]);
        hold off;
        axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])                 
        drawnow
    end
%     if visualization == 1
%         rect_position_vis = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
%         im_to_show = double(im)/255;
%         if size(im_to_show,3) == 1
%             im_to_show = repmat(im_to_show, [1 1 3]);
%         end
%         if frame == 1
%             fig_handle = figure('Name', 'Tracking');
%             imagesc(im_to_show);
%             hold on;
%             rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
%             text(10, 10, int2str(frame), 'color', [0 1 1]);
%             hold off;
%             axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
%         else
%             resp_sz = round(sz*currentScaleFactor*scaleFactors(scale_ind));
%             xs = floor(old_pos(2)) + (1:resp_sz(2)) - floor(resp_sz(2)/2);
%             ys = floor(old_pos(1)) + (1:resp_sz(1)) - floor(resp_sz(1)/2);
%             sc_ind = floor((nScales - 1)/2) + 1;
%             
%             figure(fig_handle);
%             imagesc(im_to_show);
%             hold on;
%             resp_handle = imagesc(xs, ys, fftshift(response(:,:,sc_ind))); colormap hsv;
%             alpha(resp_handle, 0.5);
%             rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
%             text(10, 10, int2str(frame), 'color', [0 1 1]);
%             hold off;
%         end
end
[~, results] = get_sequence_results(seq);
results.posp =seq.posp;
disp(['fps: ' num2str(results.fps)])
end
function z = shrinkage(x, kappa)
    z = max( 0, x - kappa ) - max( 0, -x - kappa );
%     zz=max(0,abs(x)-kappa).*sign(x);
end
% calculate the product of the cell and vector
function res = cell_DotProduct_vec(Q,vec)
QQ           = cat(4,Q{:});
[adim,bdim,cdim,~]=size(QQ);
Vecx         = ones(adim*bdim*cdim,1)*vec;
Vecx         = reshape(Vecx,[adim,bdim,cdim,length(vec)]);
resp         = QQ.*Vecx;
res          =sum(resp,4); 
end
function [hnum,hden]=calTensorProduct(Dict,h)
    Dict_mat =cat(4,Dict{:});
    [adim,bdim,cdim,len]=size(Dict_mat);
    Dict_new = reshape(Dict_mat,[adim*bdim*cdim len]);
    hnum     = Dict_new'*h(:);
    hden     = Dict_new'*Dict_new;
end