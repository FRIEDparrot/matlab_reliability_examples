%% %%%%%%%%%%%%%% 求解以下功能函数的全局灵敏度, 其中每个变量相互独立且服从(-\pi, \pi)上的均匀分布
clear, clc;
n = 3;
x_min =  -pi * ones(1, n); x_max = pi * ones(1,n);   % 定义均匀分布的最值选取
g = @(x) sin(x(:,1)) + 5 .* sin(x(:,2)).^2 + 0.1 .* x(:,3).^4 .* sin(x(:,1));

num_MCS = 1000;
% xp  = qrandstream('sobol',n,'Skip', 1e3, 'Leap',1e2); % xp(i,:) = unifinv
% (x(:,i), x_min(i), x_max(i))

% 估计失效概率

num_GRE = 100;      % Global Reliability Evaluation

%% %%%%%%%%%%%%%%%%%%%%%%%% 全局灵敏度分析 %%%%%%%%%%%%%%%%%%%%
xp = rand(num_MCS ,n) .* (x_max - x_min) + x_min;   % histogram(xp(:,3),20);
Pf = mean(g(xp) < 0);

% 在全局灵敏度分析中， 直接利用MCS分析中的样本进行修改列后进行分析即可
sensitivity_g = zeros(1, n);
for i = 1:n
    x_i = rand(num_GRE,1).* (x_max(i) - x_min(i)) + x_min(i);  % 在取值范围内进行随机抽样
    Pf_i = zeros(num_GRE,1);

    xp_pre =  xp;  % 取原先的xp -> 由于是均匀分布函数, 因此直接修改即可 
    for j = 1:num_GRE
        xp_pre(:,i) = x_i(j);
        Pf_i(j,:) = mean(g(xp_pre) < 0);
    end
    sensitivity_g(i) = 0.5 * mean(abs(Pf_i - Pf));  %% 重点是这一句, 求解的是减去Pf后的均值而不是先求均值
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 蒙特卡洛法求解对应的灵敏度 %%%%%%%%%% 
MCS_solu([0,0,0],  , g, 1000)
