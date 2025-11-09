function [Pf, Pf_mu, Pf_sigma, msc]= IMS_result(xpt, I_f, f_x, h_x, fX_mu, fX_sigma, alpha_)
%% %%%%%% 输出函数 [xpt, I_f_sig,f_x, h_x, fX_mu, fX_sigma, alpha_] %%%%%%%%
    m = length(alpha_);
    n = size(xpt, 2);
    nums = size(xpt, 1);
    
    Pf_t = zeros(1, m);  Pf2_t = zeros(1, m); % Pf2_t备用
    Pf_mu_t = zeros(m,n); Pf2_mu_t = zeros(m,n);
    Pf_sigma_t = zeros(m,n); Pf2_sigma_t = zeros(m,n);
    
    for i = 1:m
        xp = xpt(:,:,i);               % 获取对应点xp
        fail_points = I_f(:,i);        % 失效点逻辑向量
        fail_xp = xp(fail_points,:);   % 失效点的x向量
    
        Pf_t(i) = sum(f_x(fail_xp)./h_x(fail_xp))/nums;
        Pf2_t(i) = sum(f_x(fail_xp).^2./ h_x(fail_xp).^2); % 2 表示平方, 用于计算后面项; 
    
        Pf_mu_t(i,:)    = 1/nums * sum(fX_mu(fail_points,:,i)./h_x(fail_xp), 1);  % 分量均值灵敏度
        Pf2_mu_t(i,:)   = sum(1./h_x(fail_xp).^2 .* fX_mu(fail_points,:,i).^2, 1);
        
        Pf_sigma_t(i,:) = 1/nums * sum(fX_sigma(fail_points,:,i)./h_x(fail_xp), 1); % 分量方差灵敏度
        Pf2_sigma_t(i,:)= sum(1./h_x(fail_xp).^2 .* fX_sigma(fail_points,:,i).^2, 1);
    end
    clear fail_points fail_xp
    
    Pf = sum(alpha_ .* Pf_t,2);
    msc.Pf_var = sum(1/(nums * (nums -1)) * (sum(alpha_.^2' .* Pf2_t) - nums * Pf_t.^2), 2);
    msc.Pf_cov = sqrt(msc.Pf_var)./Pf;
    % 上面两条公式正确性已经检验, 可放心使用
    
    Pf_mu     = sum(alpha_' .* Pf_mu_t ,1);
    % 公式4.4.3.3 @caution: 正确性可能有误
    msc.Pf_mu_var = 1/(nums * (nums-1)) * sum(alpha_.^2' .* Pf2_mu_t - nums * (alpha_' .* Pf_mu_t).^2 , 1);  % 注意第二项不要丢掉nums
    msc.Pf_mu_cov = sqrt(msc.Pf_mu_var)./abs(Pf_mu);
    
    Pf_sigma  = sum(alpha_' .* Pf_sigma_t, 1);
    % 公式4.4.3.3 @caution: 正确性可能有误
    msc.Pf_sigma_var = 1/(nums * (nums-1)) * sum(alpha_.^2' .* Pf2_sigma_t - nums * (alpha_' .* Pf_sigma_t).^2,1);
    msc.Pf_sigma_cov = sqrt(msc.Pf_sigma_var)./abs(Pf_sigma);
end
