clear, clc;
%% %%%%%%%%%%%%%%%% 多失效模式下相关变量的线抽样方法求解 %%%%%%%%%%%%%%% 
mu_  = [2, 3, 10];
rho_ = [ 1, 0.5 ,  0.6; 0.5, 1, 0.2; 0.6, 0.2, 1]; 
sigma_  = rho_ .* ([1, 1.5 , 2]' * [1, 1.5, 2]);      %% *** 协方差矩阵 *****
g = @(x) x(:,1) +2.* x(:,2) - x(:,3) + 5;

%% %%%%%%%%%%%%%%%%% 进行相关变量的解耦部分 %%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(mu_, 2);
% 首先变化为变量 y
% mu_y = zeros(1, n);
[A,sigma_y] = eig(sigma_);  % 求解协方差阵的特征值和特征向量
A = A';                % 使用的是 y = Ax + b -> sx = A^T sy A
sigma_yd = sqrt(diag(sigma_y)');

% --------  g(x) = g(x(u)) --------- 
mu_z = zeros(1,n);
sigma_zd = ones(1,n);
sigma_z   = diag(sigma_zd.^2);
% 以u(x)为变量, 由于 y = A^T(X - mu_x), 同时 u = (y - mu_y)/sigma_y
g_new = @(u) g((A' * (sigma_yd'.* u') + mu_')');  % equation 6.2.3.4 

Pf = LS_solu_test(mu_z, sigma_z, g_new, 4e3);

function Pf = LS_solu_test(mu_, sigma_, g, num_LS)
%LS_SOLU 获取功能函数在输入变量均为独立正态分布下的线抽样法的失效概率
    % 求解归一化的设计点
    [x_ii, ~, ~] = AFOSM_solu(mu_, sigma_, g); % obtain design point by AFOSM method
    sigma_d = sqrt(diag(sigma_))'; x_i = (x_ii - mu_)./sigma_d;  clear xii;

    % use standarlized variable to derive the  important sampling direction
    alpha_ = x_i./sqrt(sum(x_i.^2)); % use standarlized vector as direction
    
    % linear sampling 
    xpp = lhsnorm(mu_, sigma_, num_LS,"on");
    [g_new, xp] = StdVar_test(mu_, sigma_, g, xpp);

    clear g xpp
    % @ caution 点积和直接相乘区别, 点积需要使用相乘再相加求解, 直接相乘只是相加之后是模, 方向并不相同
    % 例如: xp_tau = xp - sum(xp .* alpha_, 2) .* alpha; 错误代码 xp_tau = xp - xp.* alpha_ 
    xp_tau = xp - sum(xp .* alpha_, 2) .* alpha_;      % perpendicular weight of xp
    xp_r = xp - xp_tau;                                % tangential weight of xp
    coef_alpha = sqrt(sum(xp_r.^2,2)).* sign(sum(xp_r .* alpha_, 2));  % just use one column to derive the coef of xp
    coef_interp = [0.3, 0.7, 1];                                                          % interpolation coefficient
    m = size(coef_interp, 2);
    
    cj_res = zeros(num_LS, 1);
    for i = 1:num_LS
        cj_interp = coef_alpha(i) * coef_interp; % c_j interpolation
		yp_interp = zeros(1, m);
        for j = 1:m
		    xp_temp = xp_tau(i,:) + cj_interp(j) .* alpha_;  % x_p calculated (alpha->important vector)
            yp_interp(j) = g_new(xp_temp);                      
        end
        clear xp_temp
        p = polyfit(yp_interp, cj_interp, 2);  % use 2 order-function for interpolation
        % @note: p returns the 3 coefficients returned;
        cj_res(i) = p(3);  % for yp = 0, derive the correspond cj using interpolation
        clear p  cj_interp
    end
    
    Pf = mean(normcdf(-cj_res));           % obtain the failure probablity
end


%% %%%%%%%%%%%% 独立变量归一化部分 %%%%%%%%%%%%%%%%%%
function [g_new, xp_new] = StdVar_test(mu_, sigma_, g, xp)
    sigma_d = sqrt(diag(sigma_))';
    g_new = @(z) g(z .* sigma_d + mu_);
    if nargin == 4
        xp_new = (xp - mu_)./sigma_d;
    end
end