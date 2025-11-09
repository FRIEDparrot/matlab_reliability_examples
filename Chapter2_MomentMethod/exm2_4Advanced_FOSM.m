clear,clc
%%%%%% 
% 某内压圆筒型容器材料为15MnV, 随机变量取为内径D, 内压强P, 壁厚t与屈服强度sigma_s
% 这些随机变量独立且均服从正态分布, 参数在下面给出
D_mu = 460; D_sig = 7;
P_mu = 20;  P_sig = 2.4;
t_mu = 19;  t_sig = 0.8;
sigma_s_mu = 392; sigma_s_sig = 31.4;
% ---------------------------------------------------------- 
syms s p d t ss ps ds ts
g(s,p,d,t) = s - p * d/(2 * t);
% ---------------------differential-------------------------
g_s_ = diff(g(s,p,d,t),s);
g_p_ = diff(g(s,p,d,t),p);
g_d_ = diff(g(s,p,d,t),d);
g_t_ = diff(g(s,p,d,t),t);
% ------- derivates for function g------- 
g_s_var = subs(g_s_, [s,p,d,t],[sigma_s_mu, P_mu, D_mu, t_mu]);
g_p_var = subs(g_p_, [s,p,d,t],[sigma_s_mu, P_mu, D_mu, t_mu]);
g_d_var = subs(g_d_, [s,p,d,t],[sigma_s_mu, P_mu, D_mu, t_mu]);
g_t_var = subs(g_t_, [s,p,d,t],[sigma_s_mu, P_mu, D_mu, t_mu]);

% ------- mu_g , sigma_g -------------------
mu_g = g(s,p,d,t);
sigma_g = sqrt(g_s_^2 * ss^2 + g_p_^2 * ps^2 + g_d_^2 * ds^2 + g_t_^2 * ts^2);

% use the formula for lambda 
lambda_s_ = -(diff(g(s,p,d,t), s) * sigma_s_sig)/sigma_g;
lambda_p_ = -(diff(g(s,p,d,t), p) * P_sig)/sigma_g;
lambda_d_ = -(diff(g(s,p,d,t), d) * D_sig)/sigma_g;
lambda_t_ = -(diff(g(s,p,d,t), t) * t_sig)/sigma_g;

mu_values = [sigma_s_mu, P_mu, D_mu, t_mu];
sigma_values = [sigma_s_sig, P_sig, D_sig, t_sig];

lambda_s_var = subs(lambda_s_,[s,p,d,t, ss,ps,ds,ts], ...
    [mu_values(1),mu_values(2), mu_values(3),mu_values(4), sigma_s_sig,P_sig, D_sig,t_sig]);
lambda_p_var = subs(lambda_p_,[s,p,d,t, ss,ps,ds,ts], ...
    [mu_values(1),mu_values(2), mu_values(3),mu_values(4), sigma_s_sig,P_sig, D_sig,t_sig]);
lambda_d_var = subs(lambda_d_,[s,p,d,t, ss,ps,ds,ts], ...
    [mu_values(1),mu_values(2), mu_values(3),mu_values(4), sigma_s_sig,P_sig, D_sig,t_sig]);
lambda_t_var = subs(lambda_t_,[s,p,d,t, ss,ps,ds,ts], ...
    [mu_values(1),mu_values(2), mu_values(3),mu_values(4), sigma_s_sig,P_sig, D_sig,t_sig]);
sprintf("lambda_init: %f, %f, %f, %f", lambda_s_var, lambda_p_var, lambda_d_var,lambda_t_var)
lambda_values = [lambda_s_var, lambda_p_var, lambda_d_var, lambda_t_var]; 

syms b_test
Xi = mu_values;  % 初始值

beta_fst = 0;

for epoch = 1:1000
    if epoch == 1 
        beta_chk = beta_fst; % 使用0作为res1
    else
        beta_chk = beta_sec; % 上一次结果作为迭代结果检查res1
        Xi = mu_values + sigma_values .* lambda_values * beta_sec;  % 重新计算设计点
        % 利用设计点获取新的lambda值
        lambda_s_var = subs(lambda_s_,[s,p,d,t, ss,ps,ds,ts], ...
            [Xi(1), Xi(2),Xi(3), Xi(4), sigma_s_sig,P_sig, D_sig,t_sig]);
        lambda_p_var = subs(lambda_p_,[s,p,d,t, ss,ps,ds,ts], ...
            [Xi(1), Xi(2),Xi(3), Xi(4), sigma_s_sig,P_sig, D_sig,t_sig]);
        lambda_d_var = subs(lambda_d_,[s,p,d,t, ss,ps,ds,ts], ...
            [Xi(1), Xi(2),Xi(3), Xi(4), sigma_s_sig,P_sig, D_sig,t_sig]);
        lambda_t_var = subs(lambda_t_,[s,p,d,t, ss,ps,ds,ts], ...
            [Xi(1), Xi(2),Xi(3), Xi(4), sigma_s_sig,P_sig, D_sig,t_sig]);
        lambda_values = [lambda_s_var, lambda_p_var, lambda_d_var, lambda_t_var]; % 更新lambda_values值
    end
    % 注意: 保留6位有效数字, 避免精度过高导致无法计算
    x_res = vpa(mu_values + sigma_values  .* lambda_values * b_test, 6); % 每一次都要重新计算 x_res
    beta_sec = vpa(min(solve(g(x_res(1), x_res(2), x_res(3), x_res(4)) == 0, b_test)), 6);
    
    sprintf("epoch: %d, res1: %f, res2 : %f",epoch, vpa(beta_chk,6), vpa(beta_sec,6))
    if (abs(beta_chk - beta_sec))<0.001
        break;
    end
end

sprintf("design point: %f, %f, %f, %f", Xi(1), Xi(2), Xi(3), Xi(4))









