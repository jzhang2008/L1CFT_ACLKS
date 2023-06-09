function [distance_precision, PASCAL_precision, average_center_location_error, average_overlap_rate] = ...
    compute_performance_measures(positions, ground_truth, video_path, distance_precision_threshold, PASCAL_threshold)

% [distance_precision, PASCAL_precision, average_center_location_error, average_overlap_rate] = ...
%    compute_performance_measures(positions, ground_truth, video_path, distance_precision_threshold, PASCAL_threshold)
%
% For the given tracker output positions and ground truth it computes the:
% * Distance Precision at the specified threshold (20 pixels as default if
% omitted)
% * PASCAL Precision at the specified threshold (0.5 as default if omitted)
% * Average Center Location error (CLE).
%
% The tracker positions and ground truth must be Nx4-matrices where N is
% the number of time steps in the tracking. Each row has to be on the form
% [c1, c2, s1, s2] where (c1, c2) is the center coordinate and s1 and s2 
% are the size in the first and second dimension respectively (the order of 
% x and y does not matter here).

if nargin < 4 || isempty(distance_precision_threshold)
    distance_precision_threshold = 20;
end
if nargin < 5 || isempty(PASCAL_threshold)
    PASCAL_threshold = 0.5;
end
max_threshold = 50;  %used for graphs in the paper
video_path=video_path(1:end-4);
if size(positions,1) ~= size(ground_truth,1),
%     disp('Could not calculate precisions, because the number of ground')
%     disp('truth frames does not match the number of tracked frames.')
       disp('make the number of tracked frames match the number of groundtruth frames.');
       text_file=dir([video_path '*_gt.mat']);
       if ~isempty(text_file)
           load([video_path text_file(1).name]);
           positions=positions(frameIndex,:);
       else
           n = min(size(positions,1), size(ground_truth,1));
           m = max(size(positions,1), size(ground_truth,1));
           step=ceil(m/n);
           positions=positions(1:step:end,:);
       end
%     return
end

%calculate distances to ground truth over all frames
distances = sqrt((positions(:,1) - ground_truth(:,1)).^2 + ...
    (positions(:,2) - ground_truth(:,2)).^2);
distances(isnan(distances)) = [];

%calculate distance precision
distance_precision = nnz(distances < distance_precision_threshold) / numel(distances);

%calculate average center location error (CLE)
average_center_location_error = mean(distances);
% calculate precision
    precisions = zeros(max_threshold, 1);
	for p = 1:max_threshold,
		precisions(p) = nnz(distances < p) / numel(distances);
    end
    %plot the precisions
	figure('UserData','off', 'Name',['Precisions - ' ])
	plot(precisions, 'k-', 'LineWidth',2)
	xlabel('Threshold'), ylabel('Precision')
%calculate the overlap in each dimension
overlap_height = min(positions(:,1) + positions(:,3)/2, ground_truth(:,1) + ground_truth(:,3)/2) ...
    - max(positions(:,1) - positions(:,3)/2, ground_truth(:,1) - ground_truth(:,3)/2);
overlap_width = min(positions(:,2) + positions(:,4)/2, ground_truth(:,2) + ground_truth(:,4)/2) ...
    - max(positions(:,2) - positions(:,4)/2, ground_truth(:,2) - ground_truth(:,4)/2);

% if no overlap, set to zero
overlap_height(overlap_height < 0) = 0;
overlap_width(overlap_width < 0) = 0;

% remove NaN values (should not exist any)
valid_ind = ~isnan(overlap_height) & ~isnan(overlap_width);

% calculate area
overlap_area = overlap_height(valid_ind) .* overlap_width(valid_ind);
tracked_area = positions(valid_ind,3) .* positions(valid_ind,4);
ground_truth_area = ground_truth(valid_ind,3) .* ground_truth(valid_ind,4);

% calculate PASCAL overlaps
overlaps = overlap_area ./ (tracked_area + ground_truth_area - overlap_area);
% calculate average overlap rate
average_overlap_rate=mean(overlaps);
% calculate PASCAL precision
PASCAL_precision = nnz(overlaps >= PASCAL_threshold) / numel(overlaps);
end