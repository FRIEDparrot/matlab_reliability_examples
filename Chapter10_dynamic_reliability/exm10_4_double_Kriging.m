%% %%%%%%%%%%%% 双层Kriging 代理模型求解时变可靠度分析问题 %%%%%%%%%%%%%%%%%%%%%% 
clear, clc;
mu_ = [3.5, 3.5];
sigma_d = [0.3, 0.3]; sigma_ = diag(sigma_d.^2);

g = @(x, t) x(:,1).^2 .* x(:,2) - 5 .* x(:,1) .* t  + (x(:,2) + 1) .* t.^2 - 20;    % 功能函数定义
fx = @(x) NPfunc(x, mu_, sigma_d);
tspan = (0:0.05:5)';

num_DAK1 = 1024;     % 第一层Kriging 模型训练需要的样本点个数。

%% %%%%%%%%%%%%%%%%%%%%%%%%% 
xp = lhsnorm(mu_, sigma_, num_DAK1); % 抽取初始样本 x
xp = repmat(xp,   [length(tspan), 1]); 
tp = repmat(tspan,[num_DAK1,      1]); % 将样本在时域上进行扩展, 与xp相对应
yp = g(xp, tp);
%% %%%%%%%%%%%%%%% initialize the Kriging Subordinatemodel Function %%%%%%%%%%%%%
n = size(mu_, 2) + 1;
regr = @regpoly0;         
corr = @corrgauss; 		 
theta0 = 0.01 * ones(1,n);
lob  =  1e-5 .* ones(1,n);
upb  = 100 .* ones(1,n);

% %%%%%%%%%%%%%%% 初始样本点选取部分 %%%%%%%%%%% 
xpp = lhsnorm(mu_, sigma_, 30, "on");  % 初始样本点取 30 个, 并从上面逐步加点
tpp = zeros(size(xpp, 1), 1);
ypp = g(xpp, 0);                      % 初始时均以 t = 0 为训练样本点

%% %%%%%% 将xpp和tpp防在一起, 首先训练外层的 g(x,t) 代理模型 %%%%%%%%%%%%%%%%%
for epoch1 = 1:1000
	[dmodel1, ~]  = dacefit([xpp, tpp],ypp,regr, corr, theta0, lob,upb);
	[mu_g, sigma_g] = dacepredict([xp, tp],dmodel1);
	EI = EI_Pred(mu_g, sigma_g);   % 求解每个点的EI, 并且将其中EI 最大的点加入集合。
	[Cr, idx] = max(EI, [], 1);    
	% 使用 max(EI) 作为 Cr准则, 判定收敛性
	if Cr <= 0.01 * abs(min(mu_g))  % 当最大值 < 0.01时,内层训练停止
		break;
	else % 将EI最大的样本加入模型中,更新模型
        fprintf("epoch1: %d, Cr: %f,  Cr_lim: %f\n", epoch1, Cr,0.01 *abs(min(mu_g)));
        xpp = [xpp; xp(idx,:)];
        tpp = [tpp; tp(idx,:)];
		% ypp = [ypp; g(xp(idx,:), tp(idx,:))];
        ypp = [ypp; yp(idx,:)];
        % 注意这两句一定要放在后面:
        xp(idx,:) = []; tp(idx,:) = []; yp(idx,:) = [];
	end
end

clear xpp ypp xp yp tp idx Cr mu_g sigma_g
%% %%%%%% 训练内层代理模型, 用于预测x一定时, 时间域的最小值 min_t(g(x)) %%%%%%%%%%%%%
xp = lhsnorm(mu_, sigma_, num_DAK1);            % 抽取样本 x
yp = zeros(num_DAK1,1); yp_test = zeros(num_DAK1,2);
for i = 1:num_DAK1
    x_tmp = repmat(xp(i,:), [size(tspan, 1), 1]);
	yp(i) = min(dacepredict([x_tmp, tspan], dmodel1));   % 代理模型预测建立
    yp_test(i,1) = yp(i);
    yp_test(i,2) = min(g(xp(i,:), tspan));
end

clear x_tmp

% 取 20 个为初始的最小样本点, 之后的样本点使用代理模型加入
xpp = lhsnorm(mu_, sigma_, 20, "on");
ypp = zeros(1000,1);
for i = 1:20
    ypp(i) = min(g(xpp(i,:), tspan),[], 1);
end

% 由于没有t, 重新初始化Krigign模型参数
n = size(mu_, 2);
regr = @regpoly0;   
corr = @corrgauss; 		 
theta0 = 0.01 * ones(1,n);		 
lob  =  1e-5 .* ones(1,n);
upb  = 100 .* ones(1,n);

% 不使用双层嵌套模型, 否则往往浪费大量计算资源
for epoch2 = 1:1000
    [dmodel2, ~] = dacefit(xpp,ypp, regr,corr, theta0, lob, upb);
    [mu_g, sigma_g] = dacepredict(xp, dmodel2);   % 只需给出 xp, 获取对应的最小值
	
    % 根据CA函数进行自适应加点计算 -> 公式参考 10 - 77
    EI = normcdf(abs(mu_g) ./ sqrt(sigma_g)) .* fx(xp) .* sigma_g;
    [EI_pred, idx] = max(EI);
    
    % 另外计算收敛准则, 计算C_(K, MCS)
    CK = mean(normcdf(U_LearningFunc(mu_g,sigma_g)));
    if CK >= 0.9999 % 外层训练停止准则, 实际上是 CK > 0.999 约等于是 U_x > 2
        break;
    else % 进行自适应加点
        fprintf("epoch: %d, EI: %f, CK: %f, Pf: %f\n", epoch2, EI_pred, CK, mean(mu_g < 0));
        xpp = [xpp; xp(idx,:)];
        ypp = [ypp; yp(idx,:)];
        xp(idx,:) = [];
    end
end

pf = mean(mu_g <=0)



function f_X = NPfunc(x_i, mu_, sigma_d)
    f_X = prod(normpdf(x_i,mu_, sigma_d), 2);
end

% U 学习函数
function res = U_LearningFunc(mu_g, sigma_g)
    res = abs(mu_g)./sqrt(sigma_g);
end

% function res = CA_LearningFunc(mu_g, sigma_g, xp)
%     res = (1 - normcdf(abs(mu_g)./ sqrt(sigma_g))).* fx(xp) ;
% end


% @brief: 用于使用学习函数 E(I(t))来选取下一个更新时刻的积分函数(公式 10.75) -> 注意:第二个参数是 sigma_g, 因此必须使用 sqrt 求解方根
% @caution: 中间代入的部分是 mu_min - mu_g , 而不是 mu_g - mu_min
function res = EI_Pred(mu_g, sigma_g)
	mu_min = min(mu_g, [], 1);
	res =  (mu_min - mu_g) .* normcdf((mu_min - mu_g)./sqrt(sigma_g)) + ...
				sqrt(sigma_g) .* normpdf((mu_min - mu_g)./sqrt(sigma_g));
end

