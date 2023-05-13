function results=run_L1CFT_ACLKS(seq, res_path, bSaveImage, parameters)
% Sparse Regularized Correlation Filter Tracker with Adaptive Contextual
% Learning and Keyfiler Selection
%  Initialize path
% addpath('feature/');
% addpath('implementation/');
% addpath('utils/');
%project_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
%addpath(fullfile(project_dir, 'trackers', 'L1CFT_ACLKS'));
%initialize the position of the context patches
params.wsize    = [seq.init_rect(1,4), seq.init_rect(1,3)];
params.init_pos = [seq.init_rect(1,2), seq.init_rect(1,1)] + floor(params.wsize/2);

params.init_pos_left = params.init_pos - [seq.init_rect(1,4),0];
params.init_pos_right = params.init_pos + [seq.init_rect(1,4),0];
params.init_pos_bottom = params.init_pos - [0, seq.init_rect(1,3)];
params.init_pos_top = params.init_pos + [0, seq.init_rect(1,3)];

params.init_pos_lefttop = params.init_pos - [seq.init_rect(1,4),-1*seq.init_rect(1,3)];
params.init_pos_righttop = params.init_pos + [seq.init_rect(1,4),seq.init_rect(1,3)];
params.init_pos_leftbottom = params.init_pos - [seq.init_rect(1,4), seq.init_rect(1,3)];
params.init_pos_rightbottom = params.init_pos + [-1*seq.init_rect(1,4), seq.init_rect(1,3)];

params.width = seq.init_rect(1,4);
params.height = seq.init_rect(1,3);

params.keystep = 8;
params.Period = 2;
params.local_nums=4;
params.weight_const = 0.5;       %不同的取值可以决定是否考虑干扰,如果设为0，表示不考虑干扰
params.adaptive_select =1;

w = params.width;
h = params.height;
L = [w h];
[L_min,~] = min(L);
cross_l = sqrt(w*w+h*h);

yta = 0.28;
params.yta1 = yta * L_min / w;
params.yta2 = yta * L_min / h;
params.yta3 = yta * L_min / cross_l;
% HOG feature parameters
hog_params.cell_size=4;
hog_params.nDim = 31;
% Grayscale feature parameters  
grayscale_params.colorspace='gray';
grayscale_params.cell_size = 4;
grayscale_params.nDim = 1;
grayscale_params.useForColor=true;
grayscale_params.useForGray =true;
%intensity feature parameters
ic_params.tablename = 'intensityChannelNorm6';
ic_params.useForColor = false;
ic_params.cell_size = 4;
% Color name feature papameters
cn_params.tablename = 'CNnorm';
cn_params.useForGray = false;
cn_params.cell_size = 4;
cn_params.nDim = 10;

% temp = load('w2crs');
% colorname_params.w2c = temp.w2crs;
% colorname_params.nDim =10;
% colorname_params.useForColor=true;
% colorname_params.useForGray =false;
% Global feature parameters 
params.t_features = {
    struct('getFeature',@get_table_feature,'fparams',cn_params),...
    struct('getFeature',@get_table_feature,'fparams',ic_params),...
    struct('getFeature',@get_colorspace, 'fparams',grayscale_params),...  % Grayscale is not used as default
    struct('getFeature',@get_fhog,'fparams',hog_params),...
};
% Global feature parameters
params.t_global.cell_size = 4;                  % Feature cell size
params.t_global.cell_selection_thresh = 0.75^2; % Threshold for reducing the cell size in low-resolution cases
% params.t_global.normalize_power = 2;    % Lp normalization with this p
% params.t_global.normalize_size = true;  % Also normalize with respect to the spatial size of the feature
% params.t_global.normalize_dim = true;   % Also normalize with respect to the dimensionality of the feature

% Filter parameters
params.search_area_shape = 'square';    % the shape of the training/detection window: 'proportional', 'square' or 'fix_padding'
params.search_area_scale = 5;%4.5;  % 4.9;%      % the size of the training/detection area proportional to the target size
params.filter_max_area = 50^2;          % the size of the training/detection area in feature grid cells
params.min_image_sample_size = 150^2;%200^2;   % Minimum area of image samples
params.max_image_sample_size = 200^2;%250^2;   % Maximum area of image samples
params.image_sample_size = 160^2;              % Minimum area of image samples


% Learning parameters
params.padding = 1.0;         			% extra area surrounding the target
params.output_sigma_factor = 1/16;		% standard deviation for the desired translation filter output
params.scale_sigma_factor = 1/4;        % standard deviation for the desired scale filter output
params.lambda = 5*1e-2;					% regularization weight (denoted "lambda" in the paper)
params.learning_rate = 0.019;%0.95;%0.05;%0.013;%0.025;	% tracking model learning rate (denoted "eta" in the paper) 0.8是一个比较好的结果
params.number_of_scales = 5;%7;5            % number of scale levels (denoted "S"=7 in the paper)
params.scale_step = 1.01;%1.02;               % Scale increment factor (denoted "a" in the paper) 1.01
params.scale_model_max_area = 512;      % the maximum size of scale examples
params.init_strategy = 'indep';         % strategy for initializing the filter: 'const_reg' or 'indep'
params.num_GS_iter = 4;                 % number of Gauss-Seidel iterations in the learning
params.fs_rate     = 1.25;
params.distractor_thre = 0.1;           %0.15;%0.1;select the effective distractors by this parameter

params.num_scales = 33;
params.hog_scale_cell_size = 4;
params.learning_rate_scale = 0.027;%0.025;
params.scale_sigma_factor = 0.51;%1/2;
params.scale_model_factor = 1.0;
params.scale_model_max_area = 32*16;
params.scale_lambda = 1e-4;
% smoothing index
params.lr_translation = 0.88;

%ADMM parameters
params.penalty_scale_step = 5;
params.max_iterations = 2;
params.admm_iterations = 2;

% Regularization window parameters
params.reg_window_power = 10;%3.0;            % the degree of the polynomial to use (e.g. 2 is a quadratic window)
params.use_reg_window   = 1;
params.yta              = 0;

% contextual learing interval
params.interval = 1;

% Detection parameters
params.refinement_iterations = 1;       % number of iterations used to refine the resulting position in a frame
params.interpolate_response = 4;        % correlation score interpolation strategy: 0 - off, 1 - feature grid, 2 - pixel grid, 4 - Newton's method
params.newton_iterations = 5;           % number of Newton's iteration to maximize the detection scores

% Debug and visualization
params.visualization = 1;
params.debug = 0;
% params.videoPlayer  =vision.VideoPlayer;
%----------------------------------------------------------------
% params.wsize = [seq.init_rect(1,4), seq.init_rect(1,3)];
% params.init_pos = [seq.init_rect(1,2), seq.init_rect(1,1)] + floor(params.wsize/2);
params.img_files = seq.s_frames;
% Initialize
params.seq = seq;

% results= L1CFT_ACLKS(params);
results= L1CFT_ACLKS_HC(params);
% release(params.videoPlayer);
%return results to benchmark, in a workspace variable
% results.type = 'rect';
% results.res = rects;
% results.fps = fps;
% disp(['fps: ' num2str(fps)])
end