function [seq,  ground_truth] = load_video_info(video_path)

% [img_files, pos, target_sz, ground_truth, video_path] = load_video_info(video_path)

% text_files = dir([video_path '*_gt.txt']);
text_files = dir([video_path '/*_rect.txt']);
assert(~isempty(text_files), 'No initial position and ground truth (*_gt.txt) to load.')

% f = fopen([video_path text_files(1).name]);
% ground_truth = textscan(f, '%f,%f,%f,%f');  %[x, y, width, height]
% ground_truth = cat(2, ground_truth{:});
% fclose(f);
ground_truth =dlmread([video_path text_files(1).name]);
seq.format = 'otb';
seq.len = size(ground_truth, 1);
seq.init_rect = ground_truth(1,:);
%set initial position and size
% target_sz = [ground_truth(1,4), ground_truth(1,3)];
% pos = [ground_truth(1,2), ground_truth(1,1)];
%interpolate missing annotations, and store positions instead of boxes
try
%    ground_truth = [ground_truth(:,[2,1]) + (ground_truth(:,[4,3]) - 1) / 2 , ground_truth(:,[4,3])];
catch
    ground_truth=[];
end
    %list all frames. first, try MILTrack's format, where the initial and
	%final frame numbers are stored in a text file. if it doesn't work,
	%try to load all png/jpg files in the folder.
 text_files = dir([video_path '*_frames.txt']);
 if ~isempty(text_files)
     f = fopen([video_path text_files(1).name]);
     frames = textscan(f, '%f,%f');
     fclose(f);
     %see if they are in the 'imgs' subfolder or not
    if exist([video_path num2str(frames{1}, 'imgs/img%05i.png')], 'file')
        video_path = [video_path 'imgs/'];
         img_files = num2str((frames{1} : frames{2})', 'img%05i.png');
    elseif exist([video_path num2str(frames{1}, 'imgs/img%05i.jpg')], 'file')
        video_path = [video_path 'imgs/'];
         img_files = num2str((frames{1} : frames{2})', 'img%05i.jpg');
    elseif exist([video_path num2str(frames{1}, 'imgs/img%05i.bmp')], 'file')
        video_path = [video_path 'imgs/'];
         img_files = num2str((frames{1} : frames{2})', 'img%05i.bmp');
    else
         error('No image files to load.')
    end
     %list the files
     img_files = strcat(video_path,img_files);
     img_files = cellstr(img_files);
 else
     %no text file, just list all images
        video_path=[video_path 'img/'];
		img_files = dir([video_path '*.png']);
		if isempty(img_files),
			img_files = dir([video_path '*.jpg']);
			assert(~isempty(img_files), 'No image files to load.')
		end
		img_files = sort_nat({img_files.name});
        img_files = strcat(video_path,img_files);
 end
 seq.s_frames = img_files;
end

