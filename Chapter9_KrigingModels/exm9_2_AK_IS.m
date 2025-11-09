%% %%%%%%%%%%%%%%% AK-IS Kriging Model Building %%%%%%%%%%%%%%%%%%% 
% 屋架结构的可靠性分析 
clear, clc;

mu_ =     [2.9e7, 500, 1000];
sigma_d = [1.45e6, 100, 100];
sigma_ = diag(sigma_d.^2);

w = 2.4484;
t = 3.8884;
g = @(x)2.2 - 4.*100.^3./(x(:,1).* w .* t) .* sqrt((x(:,2)./w.^2).^2 + (x(:,3)./t.^2).^2);

num_AKIS_fst = 100;    % 使用50个样本点作为初始样本点
num_AKIS = 1e5;

%% %%%%%%%%%%% Kriging Model Building %%%%%%%%%%%%%%%%%%%%
[x_i, ~, ~]= AFOSM_solu(mu_, sigma_, g);
n = size(mu_,2);

%% %%%%%%%%%%% 求解重要抽样的中心抽样的概率密度 %%%%%%%%%%%%%%%%
fx_pdf = @(x)  joint_pdf(x, mu_, sigma_d);  % f_x 分布函数
hx_pdf = @(x) joint_pdf(x, x_i,sigma_d);    % 注意 h_x 的分布函数是关于x_i为中心的函数
% @note: 说明: 之后算式中使用的是 f(x)./h(x) 作为整体的函数值

xpp = lhsnorm(mu_, sigma_, num_AKIS_fst);
ypp = g(xpp);

lob = 1e-5 * ones(1,n);
upb = 20 * ones(1,n);
theta0 = 0.01 * ones(1, n);
% dmodel = dacefit(xpp, ypp, @regpoly0, @corrgauss, theta0, lob, upb);
% 这个在循环中已经调用了

%% %%%%%%%%%%% 仍然需要使用AFOSM方法构建设计点过程中的点或者使用初始选取点作为初始点构建 %%%%%%%%%%%%
xp = lhsnorm(x_i, sigma_, num_AKIS);   % 注意: 这些都是从 x_i 为中心重要抽样得到的点
% @note: 如果前面使用MCS方法后面使用IS方法, 会导致失效概率偏低, 因此预测时必须使用新获取的抽样点进行预测

cur_xp = []; cur_yp = [];
for epoch  = 1:1000
    dmodel = dacefit([xpp;cur_xp], [ypp;cur_yp], @regpoly0, @corrgauss, theta0, lob, upb);
    [mu_g, sigma_g] = predictor(xp, dmodel);  % 只是求解 xp 的预测值来计算参数
    
    U_x = U_LearningFunc(mu_g, sigma_g);
    [U_min, index_min] = min(U_x);
    if (U_min < 2)
        % 加点, 取 xp 中 U_x 最小的点加入
        sprintf("epoch: %d,min U_x: %d", epoch, U_min)
        cur_xp = [cur_xp; xp(index_min,:)];
        cur_yp = [cur_yp; g(xp(index_min,:))];
        xp(index_min,:) = [];
    else
        break;
    end 
end

%% %%%%%%%%%%%%%% 使用重要抽样计算灵敏度 %%%%%%%%%%%%%%%%%%%%%%
num_IMS = 2e4;

% 里面自带了AFOSM 求设计点的过程
% [Pf, Pf_mu, Pf_sigma,msc] = IMS_solu(mu_,sigma_,g,num_IMS);

xp_test = lhsnorm(x_i, sigma_, num_IMS);    % 用于计算失效概率
If = (predictor(xp_test, dmodel) < 0);  % 使用模型预测
Pf = sum(If.*fx_pdf(xp_test)./ hx_pdf(xp_test))  ./ num_IMS  % 失效概率
Pf_var = 1./(num_IMS -1) .* (1./num_IMS .* sum(If.* (fx_pdf(xp_test).^2 ./ hx_pdf(xp_test).^2), 1) - Pf.^2);
Pf_cov = sqrt(Pf_var)./Pf;
% 与上面的解一致即可

fX_mu = ((xp_test - mu_)./ sigma_d.^2) .* fx_pdf(xp_test);
fX_sigma = (((xp_test - mu_)./ sigma_d).^2 -1) .* fx_pdf(xp_test)./sigma_d;

Pf_mu = 1/num_IMS * sum(If .* fX_mu./hx_pdf(xp_test), 1)
Pf_sigma = 1/num_IMS * sum(If .* fX_sigma./hx_pdf(xp_test),1)

function U_x = U_LearningFunc(mu_g, sigma_g)
    U_x  = abs(mu_g)./sqrt(sigma_g);
end
% @brief: 求解多变量的联合概率密度函数(normpdf)
function f_x = joint_pdf(x, mu_,sigma_d)
    f_x = prod(1./sqrt(2 .* pi .* sigma_d.^2).* exp(-0.5 .* ((x - mu_)./(sigma_d)).^2), 2);
end

