%% Y型节点管的可靠性功能函数如下: 
% g(M, F, O) = 1 - 2e-5F  - (2e-4 * |M|)^1.2 - (2e-5 * |O|)^2.1
% 需要注意的是, 部分协方差不为零, 方差矩阵给出在 sigma_temp 中;
clc, clear;
mu_ = [2e3, 1e4,1e4];
sigma_ = diag([0.5e3, 0.2e4 , 0.4e4].^2); sigma_(1,3) = 5e5;sigma_(3,1) = 5e5;

% 初始功能函数部分
g = @(x)(1 - 2 .* 10^-5 .* x(:,2) - (2e-4 .* abs(x(:,1))).^1.2 - (2e-5 .* abs(x(:,3))).^2.1);

num_test = 10^7;
xp = lhsnorm(mu_, sigma_, num_test, 'on');   % 获取拉丁超立方抽样的点部分;
XI = find(g(xp) < 0);
Failure_Properity = size(XI, 1)/num_test;    % 计算失效概率

%% 函数中的变量变换方法: 按照 A 矩阵中的对应值进行变换, 
%% 例如, x1 =  0.9995 z1 + 0.0317 z3 + mu, 具体方法如下: 
[A,sigma_new] = eig(sigma_);  % 将矩阵进行正交分解;
x1 = @(z) mu_(1) + A(1,1) .* z(:,1) + A(1,2) .* z(:,2) + A(1,3).* z(:,3);
x2 = @(z) mu_(2) + A(2,1) .* z(:,1) + A(2,2) .* z(:,2) + A(2,3).* z(:,3);
x3 = @(z) mu_(3) + A(3,1) .* z(:,1) + A(3,2) .* z(:,2) + A(3,3).* z(:,3);
g_new = @(z)(1 - 2 .* 1e-5 .* x2(z) - (2e-4 .* abs(x1(z))).^1.2 - (2e-5 .* abs(x3(z))).^2.1);

% 变换之后得到新的分布的mu均为0;
mu_new = [0,0,0];  n = size(mu_new,2);
% 求导使用的 delta 部分
syms beta_res
delta_ = diag(repmat(0.00001, [1,3]));
sigma_new_arr = sqrt(diag(sigma_new))';

Z_i = mu_new; % 初始化值;
dPf_dz = zeros(1, n);  % 实际上是存储对于每一个z的导数值;
lambda_ = zeros(1, n);

beta_pre = 0;
for epoch = 1:1000
    for i = 1:n
        % 使用设计点计算对应的lambda值;
        dPf_dz(i) = (g_new(Z_i + delta_(i,:)) - g_new(Z_i))/0.00001;  % 存储函数g对每个变量的导数
    end
    lambda_  = -(dPf_dz .* sigma_new_arr) ./ sqrt(sum(dPf_dz.^2 .* sigma_new_arr.^2));
    % 使用新的lambda确定Zi的值,用于解方程;
    Z0  = mu_new + sigma_new_arr .* lambda_ .* beta_res;   %  注意!!!!!不要把g_new写成g
    beta_cur = double(vpasolve(g_new(Z0) == 0, beta_res));              % 求解对应的 beta 值
    % 获取新的设计点Z_i
    Z_i = mu_new + sigma_new_arr .* lambda_ .* beta_cur;
    sprintf("epoch: %d, beta: %f", epoch, beta_cur)
    if (epoch >=2 && abs(beta_pre - beta_cur)/(beta_pre)< 0.001)
        break;
    else
        beta_pre = beta_cur;
    end
end
Failure_Resolution = normcdf(-beta_cur);
disp(Failure_Resolution);

%%%%%  这个和书上结果不同， 书上结果是3.5984 e-4, 这个是3.9957e-4;

