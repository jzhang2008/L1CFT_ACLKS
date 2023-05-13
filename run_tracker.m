
% run_tracker.m

close all;
% clear all;
% Add paths
setup_paths();
% addpath('feature/');
% addpath('implementation/');
% addpath('utils/');
%choose the path to the videos (you'll be able to choose one with the GUI)
base_path = 'sequences/';

%ask the user for the video
video_path = choose_video(base_path);
index_dir=strfind(video_path,'/');
video = video_path(index_dir(end-1)+1:index_dir(end)-1);
if isempty(video_path), return, end  %user cancelled
[seq, ground_truth] = load_video_info(video_path);


% params.init_pos = floor(pos) + floor(target_sz/2);
% params.wsize = floor(target_sz);
% params.img_files = img_files;
% params.video_path = video_path;
% params.seq      = seq;

% [positions, fps] = dsst_L1(params);
results = run_L1CFT_ACLKS(seq);
% calculate precisions
positions = results.res;
thresholdSetOverlap = 0: 0.05 : 1;
success_num_overlap = zeros(1, numel(thresholdSetOverlap));
res = calcRectInt(ground_truth, positions);
for t = 1: length(thresholdSetOverlap)
    success_num_overlap(1, t) = sum(res > thresholdSetOverlap(t));
end
cur_AUC = mean(success_num_overlap) / size(ground_truth, 1);
OP_vid = sum(res >= 0.5) / numel(res);
FPS_vid = results.fps;
display([video  '---->' '   FPS:   ' num2str(FPS_vid)   '    op:   '   num2str(cur_AUC)]);
[distance_precision, PASCAL_precision, average_center_location_error,average_overlap_rate] = ...
    compute_performance_measures(positions, ground_truth,video_path);

fprintf('Center Location Error: %.3g pixels\n Overlap Rate: %.3g \n Distance Precision: %.3g %%\nOverlap Precision: %.3g %%\nSpeed: %.3g fps\n', ...
    average_center_location_error, average_overlap_rate,100*distance_precision, 100*PASCAL_precision, results.fps);