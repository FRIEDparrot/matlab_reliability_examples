%% 包络函数(Envelope Function Method)方法求解时变可靠性的示例代码
% 10.12  时变功能函数  X1^2  X2 - 5 * X1 * t  + (X_2 + 1) * t^2 - 20 的算例  
clear, clc;
mu_ = [3.5, 3.5];
sigma_d = [0.3, 0.3];
sigma_ = diag(sigma_d.^2);
tspan = (0:0.05:5-0.05)';
g = @(x, t) x(:,1).^2 .* x(:,2) - 5 .* x(:,1) .* t  + (x(:,2) + 1) .* t.^2 - 20;   % 功能函数定义
% dgdt = @(x,t) -5 .* x(:,1) +  2 .* t .* (x(:,2) + 1);                   % 功能函数对时间的导数值
% dgdx = @(x,t) [2.* x(:,1).* x(:,2)-5.*t, x(:,1).^2 + t.^2];             % 对于每个 x 值的偏导数
% dg_dxdt = @ (x,t)  [-5,  2.* t];                                        % 对于dxdt的二阶导数

% @bug: 下面的求解出是0.0000
% mu_ = [53, 122, 66.5,  100]; sigma_d = [0.1, 0.1, 0.1, 0.1]; sigma_ = diag(sigma_d.^2);
% tspan = deg2rad(95.5:0.1:155.5);
% D = @(x,t) -2 .* x(:,1) .* x(:,3) .*sin(t);
% E = @(x,t) 2 * x(:,3) .*((x(:,4) - x(:,1) .* cos(t)));
% F = @(x,t) x(:,2).^2 - x(:,1).^2 - x(:,3).^2 - x(:,4).^2 + 2 .* x(:,1) .* x(:,4) .* cos(t);
% g = @(x,t) deg2rad(0.8) - abs((deg2rad(76) + deg2rad(60) .* sin(3./4.* (t - deg2rad(95.5))) -2.* atan((D(x,t) + sqrt(D(x,t).^2 + E(x,t).^2 - F(x,t).^2))./(E(x,t)+ F(x,t)))));

%% %%%%%%%%%%%%%%% |包络函数方法求解部分| %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 数值方法求解导数大小 
n = size(mu_, 2);
delta_t = 1e-5; delta_x = diag(1e-5 * ones(1,n));
dgdt = @(x,t) (g(x,t + delta_t) - g(x,t))/delta_t;
dgdx = @(x,t) diffgx(x,t, g, delta_x);
dg_dxdt = @(x,t) diffgxt(x,t, g, delta_x, delta_t);

% 构造仅与时间有关的向量 b, b'
b0 = @(t) g(mu_, t);    b0_ = @(t) dgdt(mu_, t);                     % 对应b0处的导数值
bi  = @(t)  dgdx(mu_, t) .* sigma_d;    	 		                 % 10.3.2.1 公式中, 求解向量b
bi_ = @(t) dg_dxdt(mu_, t) .* sigma_d;                               % 求解b的导数向量
f = @(t) b0_(t) -  b0(t) .* (bi_(t) * bi(t)') ./  (bi(t) * bi(t)');  % 求解包络函数L , 并以其零点作为选用的点

%% %%%%%%%%%%%%%%% 求解包络函数的展开点 u*的位置,对应的f即为包络面  %%%%%%%%%%%%%%%%%%%%
% 利用公式求解 u*(t)的零点, 取其中不重复的作为包络面的设计点
x = [];
for i = 1:size(tspan,1)
    t0 = tspan(i);
    opts =  optimset('Display','off');
    [res, ~, exitflag, ~]= fzero(f,t0, opts);
    if (exitflag == 1)
        x = [x; res];
    else % failed to find solution
    end
end
clear res feval exitflag output

% 对y中的时间点进筛选, 去除求解到的重复时间点和时间范围以外的部分
x = round(x, 5);    % 保留5位小数, 方法1 roundn = @(x,n) 10.^n .* round(x/10.^n); roundn(x, -5)
x = x(x >= tspan(1) & x <= tspan(end));
x = unique(x,"rows","stable");   % 第三个参数stable保证unique之后, 顺序不会改变


% 求解协方差矩阵的 mu 和 Sigma 对应的值, 并且记为 Mu_L 和 Sigma_L
if ~isempty(x)    
    Mu_L = [b0(tspan(1));   b0(x); b0(tspan(end))];  % 将解得的bi(x) = dg_dx排成列阵
    ti = [tspan(1); x; tspan(end)];               % 对应的 t
    Mu_L = unique(Mu_L,"rows","stable"); ti = unique(ti,"rows","stable");
else
	Mu_L = [b0(tspan(1)); b0(tspan(end))];
    % 求解失效包络面对应的bi, 方法是
	ti = [tspan(1); tspan(end)];
end
Sigma_L = bi(ti) * bi(ti)';

% 计算每个包络点的瞬时失效概率
Pf_arr = zeros(size(ti,1),1);
for i = 1:size(ti,1)
    Pf_arr(i) = 1 - normcdf( b0(ti(i))./sqrt(bi(ti(i)) * bi(ti(i))'));  % 下面是b(ti)的模 
end

[Pf_arr, idx] = sort(Pf_arr, 1, "descend");  ti = ti(idx,:); % 按照从大到小的方法排 Pf, ti 
r_tmp = rank(Sigma_L); Pf_arr = Pf_arr(1:r_tmp,:); ti = ti(1:r_tmp,:); % 剔除失效概率小的点

%% %%%%%%%%%%%%%%%%% 构造多元正态分布函数, 求解失效概率 %%%%%%%%%%%%%%%%%%%%
num_EFM = 1e5;
xp = mvnrnd(Mu_L', Sigma_L, num_EFM);          % 使用多元正态概率函数生成三个点下的随机样本点
Pf = size(find(min(xp,[],2) < 0), 1)/num_EFM;  % 取任意一个失效, 则结构都失效

fprintf("Pf is %f", Pf); 

%% %%%%%%%%%%%%%%%% 包络函数展开点 %%%%%%%%%%%%%%%%%%%%%%%%%%% 
L  = @(u,t) b0  + bi  * u';           % 求解包络函数 L(x, t) = 0 用的函数部分
L_ = @(u,t) b0_ + bi_ * u';           % 求解 L'(u, t);


% @brief 用于求解dgdx
% @param x,t : x为1维行向量(变量个数n), t为 m 维列向量
% @param delta_x, delta_t: delta_x 是 n * n 对角矩阵, delta_t 是较小的数
% dgdx = @(x,t) (g(x + delta_x .*[1,0],t) - g(x,t))/delta_x;
function dgdx = diffgx(x,t, g, delta_x)
    n = size(x, 2);
    m = size(t, 1);
    dgdx = zeros(m,n);
    for i = 1:n
        dgdx(:,i) = (g(x + delta_x(i,:),t) - g(x, t))./delta_x(i,i);  % calculate 3 lines at a time 
    end
end

% @brief 求解dg_dxdt 
% @param x,t : x为1维行向量(变量个数n), t为 m 维列向量
% @param delta_x, delta_t: delta_x 是 n * n 对角矩阵, delta_t 是较小的数
% @note: 一般取 delta_x, delta_t > 10-5为宜, 较小往往容易出错
function dg_dxdt = diffgxt(x,t, g, delta_x, delta_t)
    n = size(x, 2);
    m = size(t, 1);
    dg_dxdt = zeros(m,n);
    % ref from computational fluid dynamics
    for i = 1:n
        dg_dxdt(:,i) = (g(x + delta_x(i,:), t+delta_t) + g(x - delta_x(i,:),t-delta_t)...
            - g(x + delta_x(i,:), t-delta_t) - g(x - delta_x(i,:), t+delta_t))./(4 * delta_x(i,i) *delta_t);
    end
end
