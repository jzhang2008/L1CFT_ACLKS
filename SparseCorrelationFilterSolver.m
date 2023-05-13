function [WW]=SparseCorrelationFilterSolver(hf_den,hf_num,lambda,rol,Iter)
%objectfunction
% min_w 1/2||Xw-y||^2+lambda||w||_1
%usage
%[]=SparseCorrelationFilterSolver(hf_den,hf_num,lambda,rol)
%INPUTS
%
%
%output:
%
%
%copyrighted by jzhang
% email:jizhangjian@sxu.edu.cn
if(nargin<5),Iter=100;end
if(nargin<4),rol=10;end
if(nargin<3),lambda=0.01;end
[row,col,dim]=size(hf_num);
beta=zeros(row,col,dim);
t   =zeros(row,col,dim);
res =zeros(1,Iter);
RELTOL  = 0.5*1*1e-2;
lambda  = repmat(lambda,[1 1 dim]);
for k=1:Iter
    f_beta=fft2(beta);
    f_t   =fft2(t);
    f_w   =bsxfun(@times,(hf_num-f_beta+rol*f_t),1./(hf_den+rol));
    t_new =shrinkage(real(ifft2(f_w)+ifft2(f_beta)/rol),lambda/rol);
%     t_new =L12(real(ifft2(f_w)+ifft2(f_beta)/rol),lambda/rol);
    beta  =beta+rol*(ifft2(f_w)-t_new);
    if max(max(sum( abs(ifft2(f_w)-t_new),2)))<RELTOL||k>=Iter  %&& max(max(sum( abs(ifft2(f_w-f_w_old)),2)))<RELTOL)
        break;
    end
    res(k)=max(max(sum( abs(ifft2(f_w)-t_new),2)));
    t=t_new;
    f_w_old=f_w;
end
WW=fft2(t_new);
end
function z = shrinkage(x, kappa)
    z = max( 0, x - kappa ) - max( 0, -x - kappa );
%     zz=max(0,abs(x)-kappa).*sign(x);
end
function z=L12(F,lambda)
[row,col,dim]=size(F);
F            =reshape(F,row*col,dim);
lambda       =reshape(lambda(:,:,1),row*col,1);
z=max(0,1-lambda./sqrt(sum(F.*F,2)))*ones(1,dim).*F;
z=reshape(z,row,col,dim);
end

