clear, clc;
% 下面为失效函数的定义, 定义为并联失效模式的情况，即取系统的功能函数为max(g1, g2); 
% 有两个变量, 其系数已经在下面给出, 注意有相关系数0.4;
g1 = @(x) 2 - x(:,2) + exp(-0.1 .* x(:,1).^2) + (0.2.* x(:,1)).^4;
g2 = @(x) 4.5 - x(:,1) .* x(:,2);
g = @(x) max(g1(x), g2(x));

mu_ = [0.85, 0];
sigma_ = [(1/3)^2, 0.4 * 1/3; 0.4 * 1/3, 1^2];   % 此处sigma_为协方差矩阵;

[Pf,Pf_mu_X, Pf_D_X] = MCS_Indirect_solu(mu_,sigma_,g);

sprintf("fail properties is %f",vpa(Pf,8))
disp("mean sensitivity:")
disp(Pf_mu_X);
disp("variance sensitivity:")
disp(Pf_D_X);

%% 说明: 程序结果与算例存在一定误差, 可能需要进一步调试