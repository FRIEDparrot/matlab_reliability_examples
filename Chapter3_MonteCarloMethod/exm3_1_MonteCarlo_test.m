% 考虑功能函数 g(X) = 4X1 -3.9998X2  + 4X3 - X4
% 均值向量和方差向量在下面给出定义
clear, clc;
g = @(x) 4 .* x(:,1) - 3.9998 .* x(:,2) + 4.*x(:,3) - x(:,4);
mu_ = [83.5, 83.5, 83.5, 150];
cov_ = [0.12, 0.12 , 0.12, 0.25];
delta_x = 0.0001;
sigma_d = cov_ .* mu_;
sigma_ = diag(sigma_d.^2);

% 一次二阶矩法(FOSM)和AFOSM的精确解为:Pf = 0.009844, 
% mean sensitivity: -0.001333, 0.001333, -0.001333, 0.000333
% variance sensitivity: 0.001579, 0.001579, 0.001579, 0.000369 

%%%%%%%% 下面可以使用  [a, b,c] = MCS_solu(mu_,sigma_, g,1e7) 来代替 %%%%%%%%%%%%%
%% *************** 蒙特卡洛方法求解部分 **************
n = size(mu_, 2);

num_MCS = 1e7;
xp = lhsnorm(mu_, diag(sigma_d.^2), num_MCS, "on");  % 拉丁超立方抽样

fail_points = find(g(xp) < 0);
Pf = size(fail_points, 1)/num_MCS;
Pf_Sigma = sqrt((Pf - Pf^2)/(num_MCS -1));
Pf_Cov = Pf - Pf_Sigma;
sprintf("failure property: %f, Cov:%f", Pf, Pf_Cov)

%% **** 求解均值可靠性局部灵敏度和方差可靠性灵敏度 **** %% 
% 首先取出抽样点中的全部失效点, 使用独立正态分布变量公式计算均值和方差
sample_ = xp(fail_points, :);  % 失效点数据

mean_sensitivity = zeros(1,n);
variance_sensitivity = zeros(1, n);
for i = 1:4
    temp = zeros(1,n); temp(i) = delta_x;
    mean_sensitivity(i) = sum((sample_(:,i) - mu_(i))./(sigma_d(i).^2))/num_MCS;  % 计算 dg_dx_(i) , 对应公式3.3.3
    variance_sensitivity(i) = sum(1./sigma_d(i) .*(((sample_(:,i) - mu_(i))./sigma_d(i)).^2 -1))/num_MCS;
end

disp(vpa(mean_sensitivity, 6));
disp(vpa(variance_sensitivity, 6));
