%% %%%%%%%%%% 基于重要抽样IS方法的预测失效精度检测 %%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc;

test_size = 4;   % 使用3组前期计算的模型数据作为对比
models = cell(1, test_size);
for i = 1:test_size - 1
     load(join(["MIAK_results\MIAK_result_backup_epoch", num2str(i * 300), ".mat"],''), "dmodel");
     models{i} = dmodel;
     clear dmodel
end
load("MIAK_results\IMS_result_final1600.mat", "dmodel");
models{test_size} = dmodel;  % 最终建立的dmodel
clear i dmodel 

%% %%%%%%%%%% 通过重要抽样方法计算不同抽样次数情况下预测的失效概率 %%%%%%%%%%
% 求解重要抽样密度分布函数 hx_pdf 部分
mu_     = [2.5, 2.5, 2.5, 2.5, 2.5];
sigma_d = [0.5, 0.5, 0.5, 0.5, 0.5];
sigma_  = diag(sigma_d.^2);
n = size(mu_, 2);

load("E:\workpack\Matlab\MATLAB_reliability_engineering\Graduate_Projects\Meta_IS_AK_SensAnalysis\MIAK_results\MIAK_result_backup_epoch20.mat", "dmodel");
Pf_epsilon = mean(fx_pred(lhsnorm(mu_, sigma_, 1e5), dmodel));
fx_pdf = @(x) joint_pdf(x,mu_, sigma_d);
px_pdf = @(x) fx_pred(x, dmodel);   % 其中, 设修正项为 pi(x); 
hx_pdf = @(x) px_pdf(x).* fx_pdf(x) ./ Pf_epsilon; 
clear dmodel  px_pdf;
%% %%%%%%%%%% 使用 MCMC 方法在失效域附近进行重要抽样, 获取用于失效概率和局部灵敏度计算的点 %%%% 
x_init = [2.0, 2.0, 2.0, 2.0, 2.0];
fprintf("preparing sampling new points for calculation.....\n");
x_pred = slicesample(x_init, 800, "burnin", 1000, "thin", 5, "pdf", hx_pdf);
fprintf("sampling points finished\n");

%% %%%%%%%%%%%%%%%%%% 这一段加上检验样本的有限元结果计算 %%%%%%%%%%%%%%%%%%%%%%% 
para_names = ["DS_FLAP_HEIGHT", "DS_LO_RIB_WIDTH1", "DS_LO_RIB_WIDTH2", ...
              "DS_LO_RIB_WIDTH3","DS_LA_RIB_WIDTH"];
g = @(x) 2.5e8 - max_stress(x, para_names);
yp_real = g(x_pred);

save("MIAK_results/test_data_backup.mat","xp","yp_real");

% %%%%%%%%%%% 绘制不同样本点下不同变量的失效概率和灵敏度情况, 以及真实值的图像 %%%% 
Pf_arr = zeros(test_size + 1, 1);
Pf_mu_arr = zeros(test_size+1, n);
Pf_sigma_arr = zeros(test_size+1, n);

fprintf("calculate loss probablity...... \n");
for i = 1:test_size + 1
    if i <= test_size 
        yp = dacepredict(x_pred, models{i});
        % 分析各个模式的失效概率和灵敏度
        If = (yp <= 0);
        nums = size(x_pred,1);
        fX_mu = ((x_pred - mu_)./ sigma_d.^2) .* fx_pdf(x_pred);
	    fX_sigma = (((x_pred - mu_)./ sigma_d).^2 -1) .* fx_pdf(x_pred)./sigma_d;
        
        Pf_arr(i) = sum(If .* fx_pdf(x_pred)./ hx_pdf(x_pred), 1) ./ nums;
        Pf_mu_arr(i,:)    = sum(If .*    fX_mu./ hx_pdf(x_pred), 1) ./ nums;
        Pf_sigma_arr(i,:) = sum(If .* fX_sigma./ hx_pdf(x_pred), 1) ./ nums;
    else % 计算真值的失效概率和灵敏度
        If = (yp_real <= 0);
        fX_mu = ((x_pred - mu_)./ sigma_d.^2) .* fx_pdf(x_pred);
	    fX_sigma = (((x_pred - mu_)./ sigma_d).^2 -1) .* fx_pdf(x_pred)./sigma_d;
        Pf_arr(i) = 0.064;
        Pf_mu_arr(i,:) = [-0.15   -0.068   -0.15   -0.1   -0.005];
        Pf_sigma_arr(i,:) = [0.0834 , 0.0295 ,0.1137 , 0.0528 , 0.0096];
    end 
end

%% %%%%%%%%%%%%%%%%%%% 绘制相应的图像 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure("Name", "Loss probablity and local sensitivity analysis");
colors = ['g', 'b', 'c', 'm', 'r'];

subplot(3,1,1);
hold on
labels = ["Pf epoch 300", "Pf epoch 600", "Pf epoch 900", "Pf final model",  "Pf real"];

for i = 1:test_size+1    
    x = labels(i);
    bar(x, Pf_arr(i), 0.4,colors(i));
    text(i-0.1, 0.03, num2str(Pf_arr(i)));
end

subplot(3,1,2);  % 局部灵敏度计算
bar(Pf_mu_arr')
legend(["ep300", "ep600", "ep900", "final",  "real"]);
title("mean value sensitivity");

subplot(3,1,3);  % 局部灵敏度计算
bar(Pf_sigma_arr')
legend(["ep300", "ep600", "ep900", "final",  "real"]);
title("square error sensitivity");

% @brief: pi(x) 值的计算, 即返回各个项失效的概率
function px = fx_pred(xp, dmodel)
    [mu_g,sigma_g] = dacepredict(xp, dmodel); % 单个输入变量, 获取方差(用于提供hx_pdf作为抽样函数)
    px = normcdf(- mu_g ./sqrt(sigma_g));     % 直接通过分布函数取均值获取 pi(x) 的每一项
end

% @brief: 求解多变量的联合概率密度函数(normpdf)
function f_x = joint_pdf(x, mu_,sigma_d)
    f_x = prod(1./sqrt(2 .* pi .* sigma_d.^2).* exp(-0.5 .* ((x - mu_)./(sigma_d)).^2), 2);
end





