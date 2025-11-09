% 取线性功能函数为g(X) = 2 - (X1 + X2)/sqrt(2); 其中X1, X2 为相互独立的标准正态变量
% 使用重要抽样法进行求解对应的结果
clear, clc;
g = @(X) 2 - (X(:,1) + X(:,2))/sqrt(2);
mu_ = [0,0];  sigma_ = [1,0;0,1];
sigma_d = sqrt(diag(sigma_))';

num_IMS = 5e3;  % 重要抽样方法样本量
n = size(mu_, 2);
% P_f = MCS_solu(mu_, sigma_, g, 1e5)
[x_i, beta_, Pf_AFOSM] = AFOSM_solu(mu_, sigma_, g);  % 求解对应的设计点以及失效概率值

%% *********** 重要抽样法求解部分 ********************** 
% 在设计点周围获取重要抽样法的样本点, 注意是标准化之后的, 因此使用1作为抽取的方差
xp = lhsnorm(x_i, diag(ones(n,1)), num_IMS); 

fail_points = find(g(xp) < 0);  % 设计点周围的所有失效点;
fail_xp = xp(fail_points,:);    % 失效点的x值和y值;
fail_yp = g(fail_xp);           % 不使用

% 对于x1和x2, 由于其概率密度函数均为正态分布函数 -> 计算在不同抽样点的概率密度;
% 计算所有**失效点**的 f_x 和 h_x 值
f_x = prod(1/sqrt(2*pi) * exp(-1/2 .* fail_xp.^2), 2);      % 求解联合概率密度的值
h_x = prod(1/sqrt(2*pi) * exp(-1/2 .* (fail_xp - x_i).^2), 2); % 重要抽样密度函数的构造: hX(x)

P_f_IMS = sum(f_x ./ h_x)/num_IMS  % 根据4.1.1, 计算对应的失效概率 
P_f_var_IMS = 1/(num_IMS -1) * (1/num_IMS * sum(f_x.^2 ./ h_x.^2) - P_f_IMS^2)  % 可以取得极小的方差 (Vpf)
P_f_cov_IMS = sqrt(P_f_var_IMS)/P_f_IMS  % 变异系数的计算 (CovPf)

% 注意, 对于重要抽样的部分点, 则均值为x_i, 方差为1, 而fX_mu是使用原始均值和方差
fX_mu = f_x .* (fail_xp - mu_ ) ./(sigma_d); % 首先求解式4.2.4
Pf_mu_IMS = sum(fX_mu ./ h_x, 1) / num_IMS % 将4.2.4代入4.2.1 获取 Pf 对mu的灵敏度

% 方差灵敏度 
fX_simga = f_x ./ sigma_d .* ((fail_xp - mu_).^2 ./ (sigma_d.^2) - 1);
Pf_sigma_IMS = sum(fX_simga ./h_x, 1)/num_IMS
% 求解灵敏度方差 (均值和上面相同)
PfV_mu_IMS    = 1/(num_IMS - 1) * ( 1/num_IMS * sum((fX_mu   ./h_x).^2) - Pf_mu_IMS.^2)
PfV_sigma_IMS = 1/(num_IMS - 1) * ( 1/num_IMS * sum((fX_simga./h_x).^2) - Pf_sigma_IMS.^2)

