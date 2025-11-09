clear, clc;
mu_ = [460, 20, 19, 392];
sigma_d = [7, 2.4, 0.8, 31.4];
sigma_ = diag(sigma_d.^2);

g = @(x) x(:,4) - x(:,2).*x(:,1)./(2 * x(:,3));

[x_i, beta_, ~] = AFOSM_solu(mu_,sigma_, g);

% 使用设计点为重要抽样点
num_TIS = 3e4;
n = length(mu_);
xp = lhsnorm(x_i,sigma_,num_TIS,'on');

r2 = sum(((xp - mu_)./(sigma_d)).^2,2);         % 将其中圆内的部分进行解决P
xp(r2 < beta_^2,:) = [];   % 截断方法

fail_points = find(g(xp) < 0); % 
fail_xp   = xp(fail_points,:); % 

f_x = @(x)joint_pdf(x, mu_, sigma_d);            
h_x = @(x)Kernel_func(x, x_i, sigma_d, beta_); % h(x)是以1为方差的概率密度函数

Pf = sum(f_x(fail_xp)./h_x(fail_xp), 1)/num_TIS;
% 其他公式略去

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
