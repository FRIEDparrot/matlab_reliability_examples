clear, clc
%%%% ====== 说明 ===== 
% 机翼的九盒段结构, 由64个杆单元件42个板单元组成, 材料为铝合金, 外载荷与各个单元强度
% 均为正态随机变量且相互独立, 外载荷P的均值mu_p = 150kg,变异系数为均值与标准差比(V_p = sigma/mu) = 0.25, 
% 第 i 个单元强度 Ri 的均值和变异系数为mu_Ri = 83.5kg, V_Ri = 0.12,(i = 68,77,78) 
% 且结构的主要失效模式的极限状态函数为 
% g(R_68, R_77, R_78, P) = 4 R_68 - 3.9998 R77 + 4.0 R_78 - P 
% 采用改进的一次二阶矩方法FOSM进行可靠性分析和可靠性灵敏度分析
%%%%
mu_R = 83.5; mu_P = 150;
sigma_R = 0.12 * mu_R; sigma_P = 0.25 * mu_P;

syms a b c d sigma_a sigma_b sigma_c sigma_d
g(a,b,c,d) =  4.0 * a - 3.9998 * b +  4.0 * c - d;

mu_g = g(a,b,c,d);
sigma_g = sqrt(4.0^2*(sigma_a)^2 + 3.9998^2* (sigma_b)^2 + ...
    4.0^2 * (sigma_c)^2 + 1 * (sigma_d)^2); % 6.2253 x 103


% use the formula for lambda 
lambda_a_ = -(diff(g(a,b,c,d), a) * sigma_a)/sigma_g;
lambda_b_ = -(diff(g(a,b,c,d), b) * sigma_b)/sigma_g;
lambda_c_ = -(diff(g(a,b,c,d), c) * sigma_c)/sigma_g;
lambda_d_ = -(diff(g(a,b,c,d), d) * sigma_d)/sigma_g;


mu_values = [mu_R, mu_R, mu_R, mu_P]; % use the mean value as the initial value
% 初始时, 取均值处作为设计点
lambda_a_var = subs(lambda_a_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_values(1),mu_values(2),mu_values(3),mu_values(4), sigma_R, sigma_R, sigma_R, sigma_P]);
lambda_b_var = subs(lambda_b_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_values(1),mu_values(2),mu_values(3),mu_values(4), sigma_R, sigma_R, sigma_R, sigma_P]);
lambda_c_var = subs(lambda_c_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_values(1),mu_values(2),mu_values(3),mu_values(4), sigma_R, sigma_R, sigma_R, sigma_P]);
lambda_d_var = subs(lambda_d_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_values(1),mu_values(2),mu_values(3),mu_values(4), sigma_R, sigma_R, sigma_R, sigma_P]);

sigma_values = [sigma_R, sigma_R, sigma_R, sigma_P];
lambda_values = [lambda_a_var, lambda_b_var, lambda_c_var, lambda_d_var];

sprintf("lambda_ini: %f, %f, %f, %f", lambda_a_var, lambda_b_var, lambda_c_var,lambda_d_var)
% solve the equation that g(a,b,c,d) = 0 to reach;
% if epoch time > 1000, break

syms b_test
x_res = mu_values + sigma_values  .* lambda_values * b_test;
Xi = mu_values;

beta_fst = 0;

for epoch = 1:1000
    if epoch == 1 
        beta_chk = beta_fst; % 使用0作为res1
    else
        beta_chk = beta_sec; % 上一次结果作为迭代结果检查res1
        Xi = mu_values + sigma_values .* lambda_values * beta_sec;  % 重新计算设计点 
        % 利用设计点获取新的lambda值
        lambda_a_var = subs(lambda_a_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
             [Xi(1), Xi(2), Xi(3), Xi(4), sigma_R, sigma_R, sigma_R, sigma_P]);
        lambda_b_var = subs(lambda_b_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
             [Xi(1), Xi(2), Xi(3), Xi(4), sigma_R, sigma_R, sigma_R, sigma_P]);
        lambda_c_var = subs(lambda_c_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
             [Xi(1), Xi(2), Xi(3), Xi(4), sigma_R, sigma_R, sigma_R, sigma_P]);
        lambda_d_var = subs(lambda_d_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
             [Xi(1), Xi(2), Xi(3), Xi(4), sigma_R, sigma_R, sigma_R, sigma_P]);
    end
    lambda_values = [lambda_a_var, lambda_b_var, lambda_c_var, lambda_d_var]; % 更新lambda_values值
    x_res = mu_values + sigma_values  .* lambda_values * b_test;  % 用于解方程
    beta_sec = solve(g(x_res(1), x_res(2), x_res(3), x_res(4)) == 0, b_test);
    sprintf("epoch: %d, res1: %.10f, res2 : %.10f", epoch, beta_chk, beta_sec)
    if (abs(beta_chk - beta_sec))<0.001
        break;
    end
end

sprintf("design point: %f, %f, %f, %f", Xi(1), Xi(2), Xi(3), Xi(4))
