%%%%%%%
% 考虑 g(X) = X1 - X2 为失效函数
% X1服从于[0, 100] 区间的均匀分布, 而X2服从于参数为12.5 的指数分布
% x1 = unifpdf(x, 0, 100); x2 = exppdf(x, 12.5);
%%%%%%
clear, clc
g = @(x1, x2) x1 - x2; % g(x)功能函数
mu = [ 50, 12.5 ];    % 初始两组数据的原始均值记为mu

%%%%%%%%%%%%%%%% initial value %%%%%%%%%%%%%%
x_i = mu;   % initialize x as mu
delta_x = 0.00001;  % used for derive the value of x

syms beta__var_

% 初始化 mu 和 sigma 数组, 以及lambda_i 值
mu_fixed = zeros(1,2); sigma_fixed =  zeros(1,2); lambda_i = 0; 

Data_Epoch = zeros(1000, 7);
num_epoch = 0; 
%%%%%%%%%%%%%%%%% 进行迭代获取相应的beta值 %%%%%%%%%%%
for epoch = 1:1000
    % ----------- 首先计算等价均值和等价的标准差 ------------------------- 
    sigma_fixed(1) = normpdf(norminv(unifcdf(x_i(1), 0, 100)), 0, 1)/unifpdf(x_i(1), 0, 100);  % calculate the sigma first
    sigma_fixed(2) = normpdf(norminv( expcdf(x_i(2),   12.5)), 0, 1)/exppdf (x_i(2), 12.5);
    
    mu_fixed(1) = x_i(1) - sigma_fixed(1) * norminv(unifcdf(x_i(1),0,100));
    mu_fixed(2) = x_i(2) - sigma_fixed(2) * norminv(expcdf(x_i(2), 12.5));

    % ----------- 代入对应的均值和g = 0的方程 ---------------------------
    %%%%%%% ======== 首先将g对x求导 ========= %%%%%%%%%%%%
    dg_dx = zeros(1,2);
    dg_dx(1) = (g(x_i(1) + delta_x, x_i(2)) - g(x_i(1), x_i(2)))/delta_x;
    dg_dx(2) = (g(x_i(1), x_i(2) + delta_x) - g(x_i(1), x_i(2)))/delta_x;
    %%%%%%% ======== 计算 lambda 的值 ======= %%%%%%%%%%% 
    lambda_i = - dg_dx .* sigma_fixed ./ sqrt(sum( dg_dx.^2 .* sigma_fixed.^2));
    
    x_i_var = mu_fixed + sigma_fixed .* beta__var_ .* lambda_i; %%% 1x2 数组, 用于后面求解方程

    beta_ = double(solve(g(x_i_var(1), x_i_var(2)) == 0, beta__var_));  % 求解方程, 获取新的设计点的 beta 值
    x_i = double(mu + sigma_fixed .* beta_  .* lambda_i);    % 仅迭代一次, 用于获取设计点的坐标, double 用于更改数据类型;
    
    % ---------- 将四个数据拼接成 1 x 5 矩阵, 然后接到Data_Epoch后面 -----
    Data = [mu_fixed, sigma_fixed, lambda_i, beta_];
    Data_Epoch(epoch,: ) = Data;  % 垂直拼接
    if (epoch >= 2) && (abs(Data_Epoch(epoch,7) - Data_Epoch(epoch - 1,7)) < 0.0001)
        epoch_num = epoch;
        break;
    end
end

Data_Epoch(1:epoch_num,:)     % 显示数据部分 
