%% 首先, 考虑航空工程中的无头铆钉的铆接应力, 通过控制铆接过程中的挤压应力, 使得
% 我们将铆钉的过程看成圆柱受压之后变形填充, 并假设:
% 1. 铆接过程中, 铆钉孔不扩大
% 2. 铆钉的体积变化忽略不计
% 3. 铆接结束后, 铆钉头部为圆柱状 
% 4. 各向同性材料
%%%%% 
clear,clc;
mu_ = [5, 20, 547.2, 5.1, 5];               % d(mm), h(mm), K(MPa), D0(mm), t(mm)
sigma_d = mu_.*[0.1, 0.02, 0.01, 0.2, 0.2]; % d(mm), h(mm), K(MPa), D0(mm), t(mm)-> variance cofficient
n_SHE = 0.15;                               % 铆钉材料的硬化指数
H_ = 2.2;                                   % 墩头高度
pressure_max = 580; % sigma_sq

% 定义功能函数为sigma_sq - sigma_max;
g =  @(x) pressure_max - x(:,3) .* (log((x(:,1).^2 .* x(:,2) - x(:,4).^2 .* x(:,5))./(2 .* H_ .* x(:,1).^2))).^n_SHE;
[Pf_mcs, Pf_mu_mcs, Pf_sigma_mcs] = MCS_solu(mu_, diag(sigma_d.^2), g, 1e7)     % -> 0.0472
clear Pf_mcs Pf_mu_mcs Pf_sigma_mcs

%% %%%%%%%%%%%%%%% 计算模型部分 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
num_MCS = 3e5;           % 蒙特拉罗法的抽样点数 -> 注意: 由于使用Kriging代理模型, 可以在更少的点下达到相同精度
num_Kriging_begin = 20;  % 从30个样本点开始(初始样本点个数) -> 从200个产生complex double
n = size(mu_, 2);
Cutting_OverFlow = 100;          % 定义一个 OverFlow 用于截断学习函数过大的样本
learning_rate = 1;               % 学习率, 每次加入几个样本点

%%%%%%%%%%%%%%%%% 进行抽样, 生成对应的样本点 (MCS生成的样本点需要较多)%%%%%%%%%%%%%%%%%%%%%%%
% 注意: 这个是进行标准均匀分布的sobol抽样;
xpp = sobolset(n,'Skip',1e5); xpp = xpp(1:num_MCS,:);       % 获取(0,1)中的MCS方法的抽样点
xp = norminv(xpp, mu_, sigma_d);                            % 通过反正态函数将均匀分布映射为正态分布
%%%%%%%%%%%%%%%% 方法一、Kriging 代理模型 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 首先初始化样本集 -----> 使用最初的几个作为样本点
X = xp(1:num_Kriging_begin, :);     % 初始 Kriging 模型的样本点
xp(1:num_Kriging_begin,:) = [];   
Y = g(X);                           % 计算初始样本值(5输入, 1输出) -> 实际上是模拟有限元过程
% 如果g(X) 是已知的, 可以考虑从中剔除距离失效域较远的部分;

%%%%%%%%%% 定义 Kriging 模型中的 超参数theta 迭代的初始值, 最大值和最小值 %%%
theta_0 = 0.01 .* ones(1,n);        % 初始定义为大小为0.01, 维数n的行向量
low_bnd = 1e-5 .*ones(1,n);         % theta 最小不超过5
upp_bnd = 20 .* ones(1,n);          % theta 最大不超过20;

min_value = 0; % 初始化最小学习函数
for epoch = 1:1000   % 进行自适应学习的Kriging模型迭代, 最多加入1000个点
    % 一阶回归方程, gauss 修正模型
    [dmodel, pref] = dacefit(X, Y, @regpoly0 , @corrgauss, theta_0, low_bnd, upp_bnd);  % 使用这些进行拟合点
    [mu_g_pred, sigma_g_pred] = predictor(xp, dmodel);                     % 将其余的全部样本值放入并进行预测
    % 获得的 mu_g_pred 是预测值, sigma_g_pred 是预测方差
    % for the opts parameter, since trial sites m > 1, return m-vector with estimated MSE (min square error) 
    
    % 对初始样本进行截断, 剔除U_x过大的点, 减小计算量
    U_x = U_LearningFunc(mu_g_pred, sigma_g_pred); % 使用U函数作为自适应学习函数, 取其从小到大作为排序;
    points_delete_ = (abs(U_x) > Cutting_OverFlow);
    xp(points_delete_,:) = [];
    U_x(points_delete_,:) = [];

    if (min_value >= 2)  % 终止条件是 U >= 2
        break;
    else
        % 每次加入学习函数值低的, 加入次数取决于 learning rate
        for i = 1:learning_rate
            [min_value,min_index] = min(U_x);    % 每一次找出其中最小的元素以及下标, 每一次加入对应的index_下标
            new_points= xp(min_index,:);
            xp(min_index,: ) = [];  % 删除xpp中的对应的行
            y_temp = g(new_points);
            if (imag(y_temp)==0)
                X = [X;new_points];      % 在x中加入对应的点
                Y = [Y;y_temp]; 
            end
        end
        sprintf("epoch: %d,U_x: %f", epoch * learning_rate ,min_value)
    end
    clear mu_g_pred sigma_g_pred y_temp
end
% 当做完之后, 整体的预测模型已经出来, 显然这时的mu_g_pred和sigma_g_pred是对应的结果, find(mu_g_pred <
% 0)为失效数量;
Pf_ref = size(find(mu_g_pred<0), 1)/size(mu_g_pred,1);  % 失效概率参考值
%% %%%%%%%%%%% 使用构建好的代理模型, 重新对原始数据进行估计对应的值 %%%%%%%%%%%%%%%%%%%%%
xp_pre = norminv(xpp, mu_, sigma_d);
y_pre  = predictor(xp_pre,dmodel);     
nums = length(y_pre);

I_f = y_pre < 0;                 % 存储失效点
Pf = length(find(I_f)) / nums;   
Pf_var = (Pf - Pf^2) / (nums-1);
Pf_cov = sqrt((1-Pf) / ((nums-1) * Pf));  % 计算局部灵敏度变异系数

%% %%%%%%%%%%%% 公式 3.3.3, 3.3.4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Pf_mu = 1./nums .* sum((xp_pre(I_f,:) - mu_) ./ sigma_d.^2, 1);   % 均值灵敏度
Pf_sigma = 1./nums .* sum( 1./sigma_d .*(((xp_pre(I_f,:)-mu_) ./ sigma_d).^2 - 1));  % 方差灵敏度

function res = U_LearningFunc(mu_g_pred, sigma_g_pred)
    res = abs(mu_g_pred)./sqrt(sigma_g_pred);
end
