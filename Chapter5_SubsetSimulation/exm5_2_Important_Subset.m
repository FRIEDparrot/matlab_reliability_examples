%% %%%%%%%%% 使用重要抽样的子集模拟方法(SUS-IMS方法) %%%%%%%%%%%%%%%%%%%%%% 
clear, clc;
mu_ = [2e4, 12, 9.82e-4, 0.04, 1e11, 2e10];
sigma_d = [1400, 0.12, 5.892e-5, 0.0048, 6e9, 1.2e9];
sigma_ = diag(sigma_d.^2);

g = @(x) 0.03 - (x(:,1) .* x(:,2).^2 ./2 .* (3.81./ (x(:,4).* x(:,6))  + 1.13 ./(x(:,3) .* x(:,5))));

%% %%%%%%%%%% SUS-IMS 抽样设置参数 %%%%%%%%%%%%%%%%%%%%% 
p0 = 0.1;          % p0 is used for adjust layer numbers in sampling, had better in 0.05-0.2
num_IMS = 2.7e4;    % use lhsnorm for the first layer
% @note: p0越大, 每次抽样点数越多, 层数也越多。单层抽样点数 = p0 * num_MCMC, 越小收敛越快, 但是抽样点也少;
x_i = mu_;          % 初始抽样点
b0 = 1e5;           % 初始化b0值(较大的值)

Pf__ = []; Pf_mu__ = []; Pf_sigma__ = [];
for epoch = 1:100
    % 初始过滤, 滤除上一层 F_k-1 区域内的点;
    xp = lhsnorm(x_i, sigma_, num_IMS);  yp = g(xp);  
    If0 = (yp <= b0);  xp= xp(If0,:);  yp = yp(If0,:); nums = size(yp,1); % nums 是截断后的样本大小;
    [B,~] = sort(yp, "ascend");  b = max(B(ceil(p0 * nums),:), 0);  If = (yp <= b);
    Pr = prod(Pf__);  % 前面所有层条件概率的乘积(Prod([]) = 1);

    % 获取所有抽样点的 qx和hx, qx 是全局概率函数即fX, hx 是条件概率;
    % 其中If0部分已经去掉了, 为了防止hx = 0的出现,一般qx = 0的地方结果就是0了
    fx =  NPfunc(xp,  mu_, sigma_d); qx = fx ./Pr; hx = NPfunc(xp, x_i, sigma_d); 
    Pf_temp = sum(If.* qx./ hx, 1) / num_IMS;       %% 重要: 不是用nums, 而是使用num_IMS作为底数

    if epoch == 1
        Pr_mu = 0;
        Pr_sigma = 0;
    else
        % 注意!!! 从0-n逐行增长不能直接使用sum(由于第一行会算错)
        Pr_mu = sum(Pf_mu__ ./ Pf__, 1);   
        Pr_sigma = sum(Pf_sigma__ ./ Pf__, 1);
    end
    Pf_mu_temp = 1./num_IMS ./Pr .* sum(If./hx.*fx .*((xp - mu_)./sigma_d.^2 - Pr_mu), 1);                            %
    Pf_sigma_temp = 1./num_IMS ./Pr .* sum(If./hx.*fx .*((((xp - mu_)./sigma_d).^2 - 1)./sigma_d - Pr_sigma), 1); %
    %%%%%%%%% 获取失效域中, 概率取值最大的点来进行下一轮重要抽样 %%%%%%%%%%%%
    [~, index]  = max(If .*fx); % 求解失效样本的概率密度
    x_i = xp(index,:); b0 = b;  % 抽样完毕之后, 更新 x_i 为这次抽样中, 同时记录上一次的b值大小
	
    Pf__ = [Pf__; Pf_temp];
	Pf_mu__ = [Pf_mu__; Pf_mu_temp];
	Pf_sigma__ = [Pf_sigma__; Pf_sigma_temp];
    clear index Pf_temp Pf_mu_temp Pf_sigma_temp
    if b<=0
        break;
    end
end

Pf = prod(Pf__);
Pf_mu = sum(Pf./Pf__ .* Pf_mu__, 1);
Pf_sigma = sum(Pf./Pf__ .* Pf_sigma__,1);
clear Pf__ Pf_mu__ Pf_sigma__ 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Part %%%%%%%%%%%%%%%%%%%%%%%
% @brief: Norm Probability Function: 假设每一个变量均为标准正态概率密度函数, 获取联合概率密度函数 
% @caution: must called as f_X = @(x)probability_func(x, mu_, sigma_d);
function f_X = NPfunc(x_i, mu_, sigma_d)
    f_X = prod(normpdf(x_i,mu_, sigma_d), 2);
end

% @brief: Normal loss Function: 对于正态分布变量,求解失效概率和均值方差灵敏度
function [Pf, Pf_mu, Pf_sigma] = NLfunc(xp, If, mu_, sigma_d, nums)
	Pf       = nnz(If)/nums;             % 初始时必定为p0(一开始<b的必定是p0)
	Pf_mu    = sum(If .* (xp - mu_)./ sigma_d.^2, 1)/nums;  
	Pf_sigma = sum(If .*(((xp - mu_)./sigma_d).^2 -1), 1)./sigma_d ./nums;
end

