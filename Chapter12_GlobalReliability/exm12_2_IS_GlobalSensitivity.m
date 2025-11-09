%% %%%%%%%%%%%%%%%% 基于重要抽样和 Markov 链估计的全局灵敏度计算 %%%%%%%%%%%%%%%%% 
clear, clc;
mu_ = [460, 20, 19, 392];
sigma_d = [7, 2.4, 0.8, 31.4];
sigma_ = diag(sigma_d.^2);
g = @(x)  x(:,4) - x(:,2) .* x(:,1) ./ (2 .* x(:,3));

% 使用一次二阶矩方法求解设计点 -> 其中, 本例中使用中心点为x_i, 范围为sigma_d的函数代替重要抽样函数 h(x)
x_i  = AFOSM_solu(mu_, sigma_, g);

fx_pdf = @(x) prod(normpdf(x, mu_, sigma_d), 2);   % 提供原始抽样密度函数
hx_pdf = @(x) prod(normpdf(x, x_i, sigma_d), 2);   % 提供重要抽样密度函数

GSA1 = MCS_soluGSA(mu_, sigma_, g, 2e5)
GSA2 = IMS_soluGSA(mu_, sigma_, g, x_i, fx_pdf, hx_pdf, 4e3)


