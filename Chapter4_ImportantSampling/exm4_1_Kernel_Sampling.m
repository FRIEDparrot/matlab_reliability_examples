% @brief:使用基于核密度函数的重要抽样方法求解线性功能函数, 其中g(X) = 2 - (X1 + X2)/sqrt(2);
clear, clc, tic;
g = @(X) 2 - (X(:,1) + X(:,2))/sqrt(2);
mu_ = [0,0];  sigma_ = [1,0;0,1];
sigma_d = sqrt(diag(sigma_))';
f_X = @(x) prod(1/sqrt(2 * pi) ./sigma_d .* exp(-1/2 * ((x - mu_)./(sigma_d)).^2), 2); % 联合概率密度函数

% @brief: 定义联合概率密度 -> 即每一个带方差部分的概率密度乘积;
x0 = [1.7, 1.7];  % 在失效域中根据经验取样本点作为迭代初始点
num_MCMC = 4e3;   % 采样个数3000点
error_threshold = 0.2/(num_MCMC);  % 截断阈值参数, 如果采样到的点对应的概率密度过小就要截除;

n = length(sigma_d);
l_i  = 6 .* sigma_d * num_MCMC .^(-1/(n+4));  % 每个变量抽样对应的链长 l_i, 采取经验数值
[x_i] = Makov_Sample(x0,l_i,f_X,g,num_MCMC);  % 总抽样, use max(x_i,[], 2) for second dim

uni_num = size(unique(x_i,'rows'),1);  % 返回独立元素的个数
lambda_  = ones(num_MCMC,1);           % 自适应局部带宽因子为lamdbda(均为1)
sigma_sample_d = std(x_i);             % 求解HM样本方差;
window_width = uni_num^(-1/(n + 4));   % 窗口宽度参数

num_KIMS = 2000;                       % 重要抽样方法的抽样个数 = num_KIMS * num_KIMS2 (10000点抽样)
num_KIMS2 = 5;                         % 从每个样本中抽取num_KIMS2个值
num_KIMS_tot = num_KIMS  * num_KIMS2;
% @brief: 重要抽样法, 从失效的样本中进行均匀抽样, 每次抽样从样本点中任取一个值j
KDE_sample = zeros(num_KIMS_tot, n);
for i = 1:num_KIMS
    j = ceil(num_MCMC * rand(1));      % 任取一个失效域分布中样本点值
    % 调用核概率密度函数 Kernel_pdf 抽样
    func_temp = @(x) Kernel_pdf(x, x_i, sigma_sample_d, lambda_, window_width, j);  % 定义临时函数来表示
    KDE_sample(num_KIMS2 *(i-1) + 1:num_KIMS2 * i, :) =  slicesample(x0, num_KIMS2, pdf=func_temp);
end
clear i j func_temp

fail_points = find(g(KDE_sample)< 0); % 失效点计算
fail_xp = KDE_sample(fail_points,:);  % 所有失效点的x值
fail_yp = g(fail_xp);                 % 所有失效点的y值

% @brief: 失效概率的计算方法: 由于是使用k(x)作为重要密度抽样函数的, 为fX/kX的平均值
fX_val = zeros(length(fail_xp), 1); hX_val = zeros(length(fail_xp), 1);
for i = 1:length(fail_xp)
    xp_temp = fail_xp(i,:);
    fX_val(i) = f_X(xp_temp); % hX_val 的一般正常值在0.1-0.3左右;
    hX_val(i) = 1/num_MCMC * sum(Kernel_pdf(xp_temp, x_i, sigma_sample_d, lambda_, window_width, 1:num_MCMC) ,1);
end
clear i xp_temp
% 由于概率密度过低的点会极大影响结果大小, 因此必须剔除部分概率密度过低的点, 避免影响结果(概率截断)
error_lines = find(hX_val < error_threshold);  % 获取并剔除概率密度过低的点;
fX_val(error_lines,:) = []; hX_val(error_lines,:) = [];
fail_points(error_lines,:) = []; fail_xp(error_lines,:) = []; fail_yp(error_lines,:) = [];

%% %%%%%%%%%%%%%%%%%%%%% 结果计算部分 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pf = sum(fX_val./hX_val)/num_KIMS_tot;                      % 获取失效概率部分(公式4.1.1)
% @caution 这两项计算与例程不相同
Pf_var_KDIS = 1/(num_KIMS_tot - 1) *(sum(fX_val.^2./hX_val.^2)/num_KIMS_tot - Pf^2); % 失效概率方差 (4.1.4)
Pf_cov_KDIS = sqrt(Pf_var_KDIS)/Pf;

fX_mu = fail_xp.*fX_val;                                              % dfx_dmu 均值局部可靠性灵敏度
Pf_mu_KDIS = 1/num_KIMS_tot * sum((1./ hX_val) .* fX_mu, 1);           % 4.2.1 (结果与截断有关)
Pf_mu_var_KDIS = 1/(num_KIMS_tot-1) * (1/num_KIMS_tot .* sum(1./hX_val .* fX_mu) - Pf_mu_KDIS.^2); % 4.2.2 
Pf_mu_cov_KDIS = sqrt(Pf_mu_var_KDIS)/abs(Pf_mu_KDIS);

fX_sigma = (fail_xp.^2 - 1).* fX_val;                                 % 方差局部可靠性灵敏度部分
Pf_sigma_KDIS = 1/num_KIMS_tot * sum( 1./ hX_val .* fX_sigma, 1);     
Pf_sigma_var_KDIS = 1/(num_KIMS_tot-1) * (1/num_KIMS_tot .* sum(1./hX_val .* fX_mu) - Pf_sigma_KDIS.^2);
Pf_simga_cov_KDIS = sqrt(Pf_sigma_var_KDIS)/abs(Pf_sigma_KDIS);
toc; 

%% %%%%%%%%%%%%%%%%%%%% function Section %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% @brief: 重要抽样核概率密度函数k_x, 获取第lines行对应k_x分量的值(适用于独立样本)
function k_x = Kernel_pdf(x, x_i, sigma_sample, lambda_, windowidth, lines)
    % lines: use 1:num_MCMC for all samples
    x_iF = x_i(lines,:);  lambda_t = lambda_(lines,:);
    % 根据每一个的k_x,构建核概率密度函数, 其中k_x的每一行代表一个分量
    n = size(x_iF, 2);
    Kx = 1./sqrt(2 * pi * sigma_sample.^2) .* exp(-1./(2 .* sigma_sample.^2).*((x - x_iF)./(windowidth .* lambda_t)).^2);
    Kx = prod(Kx, 2);
    k_x = Kx./((windowidth .* lambda_t).^n);
end

% @brief: 马尔科夫链抽样, 失效集为x_i
% empirical value: l_i  = 6 .* sigma_d * num_MCMC .^(-1/(n+4)), n is num of
% variables
function x_i = Makov_Sample(x0, l_i, f_X,g,num_MCMC)
    n = length(x0);
    x_i = zeros(num_MCMC, n);  x_i(1,:) = x0;
    for i = 2:num_MCMC
        res = zeros(1, n);  % 存储抽样得到的点
        for j = 1:n
            range = [x_i(i-1,j) - l_i(j)/2, x_i(i-1, j) + l_i(j)/2];  % 马尔科夫链抽样范围
            res(j) = unifrnd(range(1),range(2),1, 1);  % 平均抽样方法获取下一个备选转移状态
        end
        if g(res) > 0  % 抽取到未失效点
            x_i(i,:) = x_i(i-1,:);   % 直接保持原状态;
            continue
        end
        % 对于其中的每一个变量, 求解转移概率
        r = f_X(res)/f_X(x_i(i-1,:)); % 决定接收或者拒绝 
        if (rand(1) < min(r))         % 用平均分布产生数据, random('Uniform', 1)
            x_i(i,:) = res;
        else
            x_i(i,:) = x_i(i-1,:);   % 直接保持原状态;
        end
    end
end

