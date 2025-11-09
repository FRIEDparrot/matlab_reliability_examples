%% %%%%%%%%%%%%%%%%%%%%  RS 加权响应面模型的有限元模型建立方案 %%%%%%%%%%%%%%%%
clear, clc;
%% %%%%%%%%%%%%%%% 6.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = @(x) exp(0.2 .* x(:,1) + 1.4) - x(:,2);
mu_ = [0,0];
sigma_d = [1,1];
sigma_ = diag(sigma_d.^2);
fprintf("\nPf:%f, %f\n\n", RS_Linear_solu(mu_, sigma_, g), RS_Quadratic_solu(mu_, sigma_, g,2));

%% %%%%%%%%%%%%%%% 6.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
mu_ = [0, 0, 0];
sigma_d = [1, 1, 1];
sigma_ = diag(sigma_d.^2);
g = @(x) -15 .* x(:,1) + x(:,2).^2 - 3.* x(:,2) +  x(:,3).^2 + 5 .* x(:,3) + 40;

fprintf("\nPf:%f, %f\n\n", RS_Linear_solu(mu_, sigma_, g), RS_Quadratic_solu(mu_, sigma_, g,2));

%% %%%%%%%%%%%%%%% 6.3 %%%%%%%%%%%%%%%%%%%%%%%%% 
mu_ = [20, 460, 19, 392];
sigma_d = [2.4, 7, 0.8, 31.4];
sigma_ = diag(sigma_d.^2);
g = @(x) x(:,4) -  x(:,1) .* x(:,2) ./(2 * x(:,3));
fprintf("\nPf:%f, %f\n\n", RS_Linear_solu(mu_, sigma_, g), RS_Quadratic_solu(mu_, sigma_, g,2));

%% %%%%%%%%%%%%%%% 6.4 %%%%%%%%%%%%%%%%%%%%%%%%% 
mu_  = [2, 3, 10];
rho_ = [ 1, 0.5 ,  0.6; 0.5, 1, 0.2; 0.6, 0.2, 1]; 
sigma_  = rho_ .* ([1, 1.5 , 2]' * [1, 1.5, 2]);      %% *** 协方差矩阵 *****
g = @(x) x(:,1) +2.* x(:,2) - x(:,3) + 5;

fprintf("\nPf:%f, %f\n\n", RS_Linear_solu(mu_, sigma_, g), RS_Quadratic_solu(mu_, sigma_, g,2));
