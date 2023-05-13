function [cf_f,g_pre_f]=train_filter(yf,x_f,xlf_context,cf_f,g_pre_f,params,seq,use_sz,small_filter_sz)
u_f = cell(1, 2);
% mumax = [0.1;10];
mumax = 10000;
for k=1:numel(x_f)
    xlf = x_f{k};
 %intialized the related parameters
  g_f        = single(zeros(size(xlf)));
  h          = zeros(size(xlf));
  h_f        = g_f;
  eta_f      = g_f;
  mu         = 1;
  Iter       = 2;
  RELTOL     = 0.5*1*1e-2;
  mu_scale_step =params.penalty_scale_step;
  %------------------------------------------------------------------------
   S_xx0 = sum(conj(xlf) .* xlf, 3);
   S_xy = bsxfun(@times, yf, conj(xlf));
if params.adaptive_select
   if params.distractor_num>0
       u_f{k}=xlf;
       for p=1:params.distractor_num
           S_xx_context{p}=sum(conj(xlf_context{p}{1,k}).* xlf_context{p}{1,k},3);
           u_f{k}  =u_f{k}+ params.weight(p)*xlf_context{p}{1,k};
       end
       S_xx = S_xx0; 
       for p=1:params.distractor_num
           S_xx = S_xx + params.weight(p)*S_xx_context{p};
       end 
       S_uu = sum(conj(u_f{k}) .* xlf, 3);
   else
      S_xx = S_xx0; 
      u_f{k}  = xlf;
      S_uu = sum(conj(u_f{k}) .* xlf, 3); 
   end
else
    %keyframe context learning
    if ((mod(seq.frame,params.keystep * params.Period) == 0)||(seq.frame==1))
        for p=1:8
            S_xx_context{p}=sum(conj(xlf_context{p}{1,k}).* xlf_context{p}{1,k},3);
        end 
        u_f{k}= xlf + params.yta1* (xlf_context{1,1}{1,k} + xlf_context{1,2}{1,k})...
               + params.yta2 * (xlf_context{1,3}{1,k} + xlf_context{1,4}{1,k}) ...
               + params.yta3* (xlf_context{1,5}{1,k} + xlf_context{1,6}{1,k} + xlf_context{1,7}{1,k} + xlf_context{1,8}{1,k}); 
        S_xx = S_xx0 + params.yta1*(S_xx_context{1} + S_xx_context{2})...
             + params.yta2 * (S_xx_context{3} + S_xx_context{4})...
             + params.yta3 * (S_xx_context{5} + S_xx_context{6} + S_xx_context{7} + S_xx_context{8});
%         u_f{k}  = xlf;
%         S_xx    = S_xx0;
        S_uu = sum(conj(u_f{k}) .* xlf, 3);
    else
       S_xx = S_xx0; 
       u_f{k}  = xlf;
       S_uu = sum(conj(u_f{k}) .* xlf, 3);
    end
end
 %-------------------------------------------------------------------------
 if (seq.frame == 1)
            g_pre_f{k} = zeros(size(xlf));
            gamma = 0;
 else
            gamma = 0.01;%10;%15;%
 end
 % pre-compute the variables
 T = prod(use_sz);
 Sg_pre_f = sum(conj(u_f{k}) .* g_pre_f{k}, 3);
 Sfx_pre_f = bsxfun(@times, u_f{k}, Sg_pre_f);
  %ADMM
  iter          = 1;
    while(iter<=params.max_iterations)
        %   solve for G- please refer to the paper for more details
           B = S_xx + (T * (mu+gamma));
           S_lx = sum(conj(u_f{k}) .* eta_f, 3);
           S_hx = sum(conj(u_f{k}) .* h_f, 3);
           g_f = (((1/(T*(mu+gamma))) * bsxfun(@times, yf, xlf)) +(gamma/(mu+gamma))*g_pre_f{k}+(mu/(mu+gamma))*h_f- (1/(mu+gamma)) * eta_f ) - ...
                    bsxfun(@rdivide,(1/(T*(mu+gamma)) * bsxfun(@times, u_f{k}, (S_uu .* yf))+gamma/(gamma+mu)*  Sfx_pre_f+mu/(gamma+mu)*bsxfun(@times, u_f{k}, S_hx)- 1/(mu+gamma) * bsxfun(@times, u_f{k}, S_lx)), B);
%            g_f = bsxfun(@times,(S_xy-T*eta_f+mu*T*h_f+gamma*g_pre_f{k}),1./(S_xx+mu*T+gamma));
%            X   = real(ifft2(g_f+ eta_f/mu));
%            h = (T/((mu*T)+ params.lambda))*real(ifft2(mu * g_f+ eta_f));
%            [sx,sy,h] = get_subwindow_no_window(h, floor(use_sz/2), small_filter_sz);
%            t = single(zeros(use_sz(1), use_sz(2), size(h,3)));
%            t(sx,sy,:) = h;
%            h_f = fft2(t);
           X = real(ifft2(mu * g_f+ eta_f));
%              if (seq.frame == 1)
                X_temp = zeros(size(X));     
                for i = 1:size(X,3)
%                  X_temp(:,:,i) = X(:,:,i) ./  (params.reg_window.^2 + mu);
                 X_temp(:,:,i)  = shrinkage(X(:,:,i),params.reg_window/mu);
                end
%                 L = 0;
%             else
%                 X_temp=X;
%                 L = max(0,1-1./(mu*numel(X)*sqrt(sum(X_temp.^2,3))));
%                 [~,b] = sort(L(:),'descend');
%                 L(b(ceil(params.fs_rate*1/params.search_area_scale^2*numel(b)):end)) = 0;
%     
%                 X_temp = repmat(L,1,1,size(X_temp,3)) .* X_temp;
%             end
           h_f = fft2(X_temp);
           %update eta_f
           eta_f=eta_f+mu*(g_f-h_f);
           % update mu
           mu = min(mu_scale_step * mu*2, mumax);  %原始�?0.1 mumax(k)0.1 10
           iter=iter+1;
    end
           % save filter
           g_pre_f{k} = g_f;
           cf_f{k}    = g_f;
%            if seq.frame==1
%              cf_f{k}    = g_f;
%            else
%              cf_f{k}    = (1-params.learning_rate)*cf_f{k}+params.learning_rate*g_f;
%            end

        %   solve for H
%            h   = shrinkage(real(ifft2(g_f)+ifft2(l_f)/mu),omega/mu);
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
        %   update mu- betha = 10.
%         mu = min(alpha * mu, mumax);
         
 end
end
 function z = shrinkage(x, kappa)
    z = max( 0, x - kappa ) - max( 0, -x - kappa );
%     zz=max(0,abs(x)-kappa).*sign(x);
end