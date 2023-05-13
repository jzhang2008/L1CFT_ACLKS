function [pos_local_trans,num_max,local_max]=generate_distractor_pos(scores_fs,small_filter_sz,trans_vec,params)
%scores_fs=double(scores_fs);
peak = max(scores_fs(:));
resp = scores_fs/peak;
resp = circshift(resp,floor([size(resp,1),size(resp,2)]/2));
% figure(3),surface(double(resp));
% colorbar
% colormap('jet')
resp = circshift(resp,round(trans_vec));
[row_res,col_res]=size(resp);
[row,col] = find(resp==max(resp(:)));
BW   = imregionalmax(resp);
% select the suitable scale factors in order to avoid beyond the boundary
BW_sz = size(BW);
scale_factor =BW_sz./small_filter_sz;
zoom_factor  = ones(size(small_filter_sz))*2;
zoom_factor(scale_factor<2.5) =1;
small_filter_sz_new =small_filter_sz.*zoom_factor;
Bys  =ceil(size(BW,1)/2-small_filter_sz_new(1)/2):ceil(size(BW,1)/2+small_filter_sz_new(1)/2);
Bxs  =ceil(size(BW,2)/2-small_filter_sz_new(2)/2):ceil(size(BW,2)/2+small_filter_sz_new(2)/2);
% Bys  =floor(size(BW,1)/2-small_filter_sz_new(1)/2):floor(size(BW,1)/2+small_filter_sz_new(1)/2);
% Bxs  =floor(size(BW,2)/2-small_filter_sz_new(2)/2):floor(size(BW,2)/2+small_filter_sz_new(2)/2);
BW(Bys,Bxs)=0;
CC = bwconncomp(BW);
local_max = zeros(1,params.local_nums);
pos_local = zeros(params.local_nums,2);
num_max   =0;
if ~isempty (CC.PixelIdxList) 
    while num_max<params.local_nums
        CC = bwconncomp(BW);
        if isempty (CC.PixelIdxList)
            break;
        end
        local_max_value = zeros(length(CC.PixelIdxList),1);
          for i = 1:length(CC.PixelIdxList)
             local_max_value(i) = resp(CC.PixelIdxList{i}(1));
          end
            lmax = max(local_max_value);
         if lmax< params.distractor_thre
             break;
         else
            num_max=num_max+1;
            [local_row,local_col ]       = find((resp.*BW)==lmax);
            pos_local(num_max,:)        =[local_row(1),local_col(1)];
%             [pos_local(num_max,1),pos_local(num_max,2)]=find(resp==lmax);
            local_max(1,num_max) = lmax;
            Bys_id = floor(pos_local(num_max,1)-small_filter_sz_new(1)/2):floor(pos_local(num_max,1)+small_filter_sz_new(1)/2);
            Bxs_id = floor(pos_local(num_max,2)-small_filter_sz_new(2)/2):floor(pos_local(num_max,2)+small_filter_sz_new(2)/2);
            %check for out-of-bounds coordinates
            Bys_id(Bys_id<1)=1;
            Bxs_id(Bxs_id<1)=1;
            Bys_id(Bys_id>row_res)=row_res;
            Bxs_id(Bxs_id>col_res)=col_res;
            BW(Bys_id,Bxs_id)=0;
        end
    end
end
pos_local_trans=pos_local(1:num_max,:)-ones(num_max,1)*[row,col];
local_max(local_max==0)=[];
%    if length(CC.PixelIdxList) > 1
%       local_max = zeros(length(CC.PixelIdxList),1);
%       for i = 1:length(CC.PixelIdxList)
%           local_max(i) = resp(CC.PixelIdxList{i}(1));
%       end
%       [local_max,local_max_idx] = sort(local_max, 'descend');
%    end
%  k_num=0;
%  
%  pos_local=zeros(params.local_nums,2);
%  while k_num<=params.local_nums
%      if local_max(1)<params.distractor_thre
%          num_max=0;
%          break;
%      end
%      k_num=k_num+1;
%      local_CC_idx = CC.PixelIdxList{local_max_idx}(1);
%  end
%    
%    if length(local_max)<params.local_nums
%       num_max = length(local_max);
%    else
%       num_max = params.local_nums;
%    end
%  
%  if num_max>0
%      pos_local=zeros(num_max,2);
%      for i=1:num_max
%          [pos_local(i,1),pos_local(i,2)]=find(resp==local_max(i));
%      end
%  end
%  pos_local_trans=pos_local-ones(num_max,1)*[row,col];