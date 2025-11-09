%% %%%%%%%%%%%%%%%%%%% 加权二次响应面法 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear, clc;
% mu_ = [0,0];
% sigma_d = [1,1];
% sigma_ = diag(sigma_d.^2);
% g = @(x) exp(0.2 .* x(:,1) + 1.4) - x(:,2);
mu_     = [0, 0, 0];
sigma_d = [1, 1, 1];
sigma_ = diag(sigma_d.^2);
g = @(x) 1./40 .* x(:,1).^4 + 2 .* x(:,2).^2 + x(:,3)+ 3;

% [Pf, Pf_mu, Pf_sigma] = MCS_solu(mu_, sigma_,g, 1e7);   % 3.09e-4  & 0.0000	-0.0000	-0.0010 
%% %%%%%%%%%%%%%%%%%%% 加权非线性响应面方法求解部分 %%%%%%%%%%%%%%%%%%% 
n = size(mu_, 2);
w = @(gx) min(abs(gx) + 1e-3)./(abs(gx) + 1e-3);   % 这个是w1
%fx = @(xp) prod(normpdf(xp, mu_, sigma_d), 2);  % 计算联合概率密度函数, 查一下源码
%w = @(xp, gx) min((abs(gx) + 1e-3)./fx) ./ ((abs(gx) + 1e-3)./fx);

f = 1;                   % 插值系数
beta_pre = 1e-3;         % 初始值,设置为 0
x_i = mu_; x0 = mu_; y0 = g(mu_); 

for epoch = 1:1000
    xp = buncher_sample(x_i, sigma_d, "f", f);  % 建立样本点
    y = g(xp);  m = length(y);            % 获取真实功能函数值 
    W = diag(w(y));                       % 权重矩阵
    A = [ones(m,1), xp, xp.^2];           
    b = ((A' * W * A) \ (A' * W  * y))';  % 最小二乘求解b 

    % 使用b获取新的预测函数值
    g_new = @(x) b(1) + sum(b(2:n + 1) .* x, 2) + sum(b(n+2: 2*n+1).* x.^2, 2);
    
    [x_i, beta_res, ~] = AFOSM_solu(mu_, sigma_, g_new);
    
    if abs((beta_res - beta_pre)/beta_pre) < 0.001
        break
    else
        x_i = mu_ + (x_i - mu_) .* (y0)./(y0 - g_new(x_i));
        sprintf("epoch: %d, beta: %f", epoch, beta_res)
        beta_pre = beta_res;
    end
end

%% %%%%%%%%%%%%%%%% 使用MCS方法, 计算响应面失效概率  %%%%%%%%%%%%%%%%%%
num_MCS = 5e6;
xpp = lhsnorm(mu_, sigma_, num_MCS,"on");
fail_points = find(g_new(xpp) < 0);
Pf = size(fail_points,1)/num_MCS;



