clear, clc; format long; 

mu_ = [2e4, 12, 9.82e-4, 0.04, 1e11, 2e10];
sigma_d = [1400, 0.12, 5.892e-5, 0.0048, 6e9, 1.2e9];
sigma_ = diag(sigma_d.^2);

g = @(x) 0.03 - (x(:,1) .* x(:,2).^2 ./2 .* (3.81./ (x(:,4).* x(:,6))  + 1.13 ./(x(:,3) .* x(:,5))));

p0 = 0.15;          % p0 is used for adjust layer numbers in sampling, had better in 0.05-0.2
% p0越大, 每次抽样点数越多, 层数也越多。单层抽样点数 = p0 * num_MCMC, 越小收敛越快, 但是抽样点也少;
num_MCS = 4.2e3;    % use lhsnorm for the first layer
num_MCMC = 30;      
f_X = @(x)NPfunc(x, mu_, sigma_d);  % 每一次以 mu 为中心进行抽样
%%%%%%%%%%%%%%%%%%%%%% Monte Carlo Solution %%%%%%%%%%%%%%%%%%%%%%%%
% [Pf, Pf_mu, Pf_sigma] = MCS_solu(mu_, sigma_, g, 1e7);
%% %%%%%%%%%%%%%%%%%%% SUS-MCMC method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(mu_, 2);
xp = lhsnorm(mu_, sigma_, num_MCS);  yp  = g(xp);

[B, ~] = sort(yp, 'ascend'); % also returns a collection of index vectors for any of the previous syntaxes
index_ = ceil(p0 * num_MCS); 
b = max(B(index_), 0);  % 查看所在位置的部分是否已经满足 < 0 的条件

If = (yp <= b);                           % 用于计算Pf_mu, Pf_sigma
Pf__       = nnz(If)/num_MCS;             % 初始时必定为p0(一开始<b的必定是p0)
Pf_mu__    = sum(If .* (xp - mu_)./ sigma_d.^2, 1)/num_MCS;
Pf_sigma__ = sum(If .*(((xp - mu_)./sigma_d).^2 -1), 1)./sigma_d ./num_MCS;

% 当不满足范围内b0小于0时,则进行分层抽样, 分层抽样中使用 MCMC方法
while b > 0 
    % @brief: 首先从范围内抽取 <= b 的部分, 即利用I取出条件概率;
    xpp = xp(If,:);   % 取出其中的 <= b 部分, 计算g值
    ypp = g(xpp);
    %% 使用MCMC方法进行继续抽样, 在每一个子集失效点下通过扩展抽样确定失效概率
    xp = zeros(size(xpp, 1) * num_MCMC, n);
	yp = zeros(size(xpp, 1) * num_MCMC,1);
    l_i = 6 .* sigma_d .* num_MCMC .^(-1/(n+4));  % copy from MCMC_Sample
    for i = 1:size(xpp,1)
        [temp_xp, temp_yp] = MCMC_Sample(xpp(i,:), l_i, f_X, g, b, num_MCMC, ypp(i)); % 以 b0 为失效准则
        xp((i-1) * num_MCMC+1: i*num_MCMC, :) = temp_xp;
        yp((i-1) * num_MCMC+1: i*num_MCMC, :) = temp_yp;
    end
    [B, ~] = sort(yp, "ascend");  length_temp  = size(yp,1);
    b = max(B(ceil(p0 * length_temp)), 0);
	
	%% %%%%%%%%%%%%  referesh the loss probability  %%%%%%%%%%%%%% 
    If = (yp <= b);
    Pf_temp = sum(If)/length_temp;
	
    % 注意: 每一项都减去, 因此要把括号包含在If中 -> 注意 sum(Pf_mu__) 之后都要加上除以\hat{P}_F, 这个是公式问题;
    Pf_mu_temp = sum(If.*( (xp - mu_)./sigma_d.^2- sum(Pf_mu__./Pf__ ,1) ) , 1)./length_temp;  % 实际上计算大小使用 /n 进行计算
	Pf_sigma_temp = sum(If .*((((xp - mu_)./ sigma_d ).^2 - 1) ./ sigma_d  - sum(Pf_sigma__./Pf__,1)) , 1)./length_temp;

    Pf__ = [Pf__; Pf_temp];
    Pf_mu__ = [Pf_mu__; Pf_mu_temp];
    Pf_sigma__ = [Pf_sigma__; Pf_sigma_temp];
    clear Pf_temp Pf_mu_temp Pf_sigma_temp xpp ypp 
    %[xp, uniq] = unique(xp, "rows"); yp = yp(uniq); % 取出其中独立的部分;
    % clear uniq temp_xp temp_yp
end

Pf = prod(Pf__, 1);
Pf_mu = sum((Pf./Pf__) .* Pf_mu__, 1);
Pf_sigma = sum((Pf./ Pf__).* Pf_sigma__, 1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Part %%%%%%%%%%%%%%%%%%%%%%%
% @brief: 假设每一个变量均为标准正态概率密度函数, 获取联合概率密度函数
% @ caution: must called as f_X = @(x)probability_func(x, mu_, sigma_d);
function f_X = NPfunc(x_i, mu_, sigma_d)
    f_X = prod(normpdf(x_i,mu_, sigma_d), 2);
end


% 优化方法: 
% 1. 由于子集模拟方法的抽样点数会随着层数增加而急剧增加, 
% 因此为了直接定义一个总的单层抽样个数和一个增益倍率(1- 1.2)左右, 可以保证维持对应的总样本池数量稳定, 提高抽样效率