function [Pf, Pf_mu, Pf_sigma,msc] = IMS_solu(mu_, sigma_, g_arr, num_MXS)
    %IMS_SOLU 重要抽样方法求解多模式下的混合概率密度(兼容单失效模式);
    %@note: 仅对独立正态变量有效, 且msc中部分项可能有误
    if (nargin == 3)
        num_MXS = 2500;
    end
    %% %%%%%%%%% 使用AFOSM方法求解设计点以及多个失效模式的可靠度 %%%%%%%
    n = length(mu_);
    m = length(g_arr);                % 失效模式个数
    sigma_d = sqrt(diag(sigma_))';    % 
    
    x_arr = zeros(m, n); beta_ = zeros(1,m);
    for i = 1:m
        [x_arr(i,:), beta_(i), ~] = AFOSM_solu(mu_, sigma_, g_arr{i});
    end
    % 定义联合概率密度 -> 即每一个带方差部分的概率密度乘积;
    %% %%%%%%%%%%%%%%%%% 多失效模式下的混合抽样 %%%%%%%%%%%%%%%%%%%%%%
    If_sig = false(num_MXS,m);       % 失效标志位
    xpt = zeros(num_MXS, n, m);       % xp_tot, 简称为xpt
    alpha_ = normcdf(-beta_)/sum(normcdf(-beta_));
    for i = 1:m
        g = g_arr{i};  % 获取对应的失效函数 -> ****不要直接使用g***
        xpt(:,:,i) = lhsnorm(x_arr(i,:), sigma_,num_MXS);  % xp, 抽样点     
        If_sig(g(xpt(:,:,i))< 0,i) = 1;                   % 取出其中失效的点
    end
    clear xp fail_xp
    %% %%%%%%%%%%%%%%%% 结果求解部分 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_x = @(x) joint_pdf(x,mu_,sigma_d);
    h_x = @(x) Kernel_func(x, x_arr, sigma_d, beta_);
    fX_mu = zeros(num_MXS,n,m);  fX_sigma = zeros(num_MXS, n,m);
    for i = 1:m
        xp = xpt(:,:,i);
        fX_mu(:,:,i) = ((xp- mu_)./ sigma_d.^2) .* f_x(xp);
        fX_sigma(:,:,i) = (((xp - mu_)./ sigma_d).^2 -1) .* f_x(xp)./sigma_d;
    end
    [Pf, Pf_mu, Pf_sigma, msc] = IMS_result(xpt, If_sig, f_x, h_x, fX_mu, fX_sigma, alpha_);
end

%% %%%%%%%%%%%%%%%%%%%%% 函数部分 %%%%%%%%%%%%%%%%%%%%%%%%%%% 
% @brief: 多变量的混合重要抽样密度函数(直接使用正态函数代替)
function h_x = Kernel_func(x, x_arr, sigma_d, beta_)
    m = length(beta_);                    % 失效模式个数
    alpha_ = normcdf(-beta_)/sum(normcdf(-beta_));  % alpha, 每个失效模式权重
    temp = 0;
    for i = 1:m
        x_i = x_arr(i,:);
        H = prod(1./sqrt(2 * pi .* sigma_d.^2).* exp(-0.5 .* ((x - x_i)./(sigma_d)).^2),2);
        temp = temp + alpha_(i) * H;
    end
    h_x = temp;
end

% @brief: 求解多变量的联合概率密度函数
function f_x = joint_pdf(x, mu_,sigma_d)
    f_x = prod(1./sqrt(2 * pi * sigma_d.^2).* exp(-0.5 .* ((x - mu_)./(sigma_d)).^2), 2);
end
