function [PSR]=CalPSR(response,reject_sz,use_sz)
resp=fftshift(response);
[row,col,~]=find(resp==max(resp(:)),1);
half_sz    =floor(reject_sz/2);
mask =ones(use_sz);
mask(row-half_sz:row+half_sz,col-half_sz:col+half_sz)=zeros(reject_sz);
residue_res=resp.*mask;
mean_val =mean(residue_res(residue_res(:)~=0));
std_val  =std(residue_res(residue_res(:)~=0));
resp_max = resp(row,col);
PSR      =(resp_max-mean_val)/std_val;


