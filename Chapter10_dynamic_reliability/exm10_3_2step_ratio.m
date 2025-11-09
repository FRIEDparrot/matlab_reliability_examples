%% %%%%%%%%%%%%%%%%%%%% 使用跨越率方法计算时变失效概率与可靠性 %%%%%%%%%%%%%%%%%%%%
% bug: 没修好AFOSM_solu问题 ->  
clear, clc;
t = (0:0.1:5)';
mu_ = [3.5, 3.5]; sigma_d = [0.3, 0.3]; sigma_ = diag(sigma_d.^2);
g = @(x, t) x(:,1).^2 .* x(:,2) - 5 .* x(:,1) .* t + (x(:,2) + 1) .* t.^2 - 20;
%% %%%%%%%%%%%%%%%%%%%% example 10.2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% t = deg2rad(95.5: 0.05*180/pi :155.5);
% mu_ = [53, 122, 66.5,  100]; sigma_d = [0.1, 0.1, 0.1, 0.1]; sigma_ = diag(sigma_d.^2);
% 
% D = @(x,t) -2 .* x(:,1) .* x(:,3) .*sin(t);
% E = @(x,t) 2 * x(:,3) .*((x(:,4) - x(:,1) .* cos(t)));
% F = @(x,t) x(:,2).^2 - x(:,1).^2 - x(:,3).^2 - x(:,4).^2 + 2 .* x(:,1) .* x(:,4) .* cos(t);
% g = @(x,t) deg2rad(0.8) - abs((deg2rad(76) + deg2rad(60) .* sin(3./4.* (t - deg2rad(95.5))) -2.* atan((D(x,t) + sqrt(D(x,t).^2 + E(x,t).^2 - F(x,t).^2))./(E(x,t)+ F(x,t)))));

% Pf = time_MCS(mu_,sigma_,g,t, 1e5);
% 先计算初始的beta0的值 -> 如果求解失败, 则这个方向上没有设计点
[x0,beta_0] = AFOSM_solu(mu_, sigma_, @(x)g(x, t(1)));  % 获取初始时的失t0,t0时刻的设计点
x0 = (x0 - mu_)./sigma_d;                               % 进行规范化变量
Pf0 = normcdf(-norm(x0));                               % 0时刻失效点在正态空间中的失效概率
delta_t = t(2) - t(1);                                  % 获取时间步进值

cross_rate = zeros(length(t)-1, 1);

for i = 1:length(t) - 1
	%先使用AFOSM方法求解 tau 和 tau + delta_tau 部分
    [x_i1, beta_1, ~, exitflag1] = AFOSM_solu(mu_, sigma_, @(x)g(x, t(i)));
    [x_i2, beta_2, ~, exitflag2] = AFOSM_solu(mu_, sigma_, @(x)g(x, t(i+1)));

    x_1 = (x_i1 - mu_) ./ sigma_d;  x_2 = (x_i2 - mu_) ./ sigma_d; % 先将设计点以mu为中心进行规范化, 获取正态空间中的设计点 x_1, x_2; (U)
	alpha_1= x_1 ./ beta_1;         alpha_2 =  x_2./beta_2;        % 分别计算 alpha 1, alpha_2 的大小, 是原点指向设计点的单位方向向量
	
    rho = - dot(alpha_1,alpha_2);  	                               % 使用两个向量点乘负值作为相关系数
    fprintf("rho:%f\n", rho);
    R = [1,rho;rho,1];               % 协方差矩阵
    if (det(R) < 5e-6)               % 奇异矩阵
        error("solution failed, try increase the time step");
    end

    % 注意: 一般的 AFOSM 求解出的 beta 会产生符号问题, 此处暂时强制纠正
	cross_rate(i) = mvncdf([beta_1, -beta_2],[0, 0], R) ./delta_t; % 跨越率 nv(tau)的求解(10-36);
end

Pf = 1 - (1-Pf0) .* exp(-sum(cross_rate).* delta_t);   % (10-34)采用小段上近似相等代替积分公式;
fprintf("Pf = %f\n",Pf)
% Pf = time_MCS(mu_, sigma_, g, t, 2e4); -> 0.1824 
