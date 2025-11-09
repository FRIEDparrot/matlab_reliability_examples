clear, clc 
%%%% ====== 说明 ===== 
% 机翼的九盒段结构, 由64个杆单元件42个板单元组成, 材料为铝合金, 外载荷与各个单元强度
% 均为正态随机变量且相互独立, 外载荷P的均值mu_p = 150kg,变异系数为均值与标准差比(V_p = sigma/mu) = 0.25, 
% 第 i 个单元强度 Ri 的均值和变异系数为mu_Ri = 83.5kg, V_Ri = 0.12,(i = 68,77,78) 
% 且结构的主要失效模式的极限状态函数为 
% g(R_68, R_77, R_78, P) = 4 R_68 - 3.9998 R77 + 4.0 R_78 - P 
% 采用一次二阶矩方法FOSM进行可靠性分析和可靠性灵敏度分析
%%%%
mu_R = 83.5; mu_P = 150;
sigma_R = 0.12 * mu_R; sigma_P = 0.25 * mu_P;

syms a b c d
g(a,b,c,d) =  4.0 * a - 3.9998 * b +  4.0 * c - d;

syms sigma_a sigma_b sigma_c sigma_d
mu_g  = g(a,b,c,d);

% 注意此处不能将sigma表示为mu的函数 %
sigma_g = sqrt(4.0^2*(sigma_a)^2 + 3.9998^2* (sigma_b)^2 + ...
    4.0^2 * (sigma_c)^2 + 1 * (sigma_d)^2); % 6.2253 x 103

beta = mu_g/sigma_g;
beta_v = vpa(subs(beta,[a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
    [mu_R,mu_R,mu_R,mu_P, sigma_R, sigma_R, sigma_R, sigma_P]));
% derive 1st derivation for g(x) ----> diff(g(a,b,c,d),a,1);

% ----- partial derivate for mu_X ----
beta_a_ = diff(beta, a);  
beta_b_ = diff(beta, b);
beta_c_ = diff(beta, c);
beta_d_ = diff(beta, d);

sprintf("beta = %f",beta_v)
sprintf("Pf = %f", normcdf(-beta_v))

Pf_beta_ = -1/sqrt(2 * pi) * exp(-mu_g^2/(2 * sigma_g^2));

sen_a = subs(beta_a_ * Pf_beta_,[a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_R,mu_R,mu_R,mu_P, sigma_R, sigma_R, sigma_R, sigma_P]);
sen_b = subs(beta_b_ * Pf_beta_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_R,mu_R,mu_R,mu_P, sigma_R, sigma_R, sigma_R, sigma_P]);
sen_c = subs(beta_c_ * Pf_beta_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_R,mu_R,mu_R,mu_P, sigma_R, sigma_R, sigma_R, sigma_P]);
sen_d = subs(beta_d_ * Pf_beta_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_R,mu_R,mu_R,mu_P, sigma_R, sigma_R, sigma_R, sigma_P]);

sprintf("mean sensitivity: %f, %f, %f, %f", sen_a, sen_b, sen_c, sen_d)

% -----
beta_sigma_a_ = diff(beta, sigma_a);  
beta_sigma_b_ = diff(beta, sigma_b);
beta_sigma_c_ = diff(beta, sigma_c);
beta_sigma_d_ = diff(beta, sigma_d);
sen_sig_a = subs(beta_sigma_a_ * Pf_beta_,[a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_R,mu_R,mu_R,mu_P, sigma_R, sigma_R, sigma_R, sigma_P]);
sen_sig_b = subs(beta_sigma_b_ * Pf_beta_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_R,mu_R,mu_R,mu_P, sigma_R, sigma_R, sigma_R, sigma_P]);
sen_sig_c = subs(beta_sigma_c_ * Pf_beta_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_R,mu_R,mu_R,mu_P, sigma_R, sigma_R, sigma_R, sigma_P]);
sen_sig_d = subs(beta_sigma_d_ * Pf_beta_, [a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
     [mu_R,mu_R,mu_R,mu_P, sigma_R, sigma_R, sigma_R, sigma_P]);

sprintf("varance sensitivity: %f, %f, %f, %f", sen_sig_a, sen_sig_b, sen_sig_c,sen_sig_d)
%sprintf()
% sprintf("beta_a_ = %f check %f", subs(beta_a_,[a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
%     [mu_R,mu_R,mu_R,mu_P, sigma_R, sigma_R, sigma_R, sigma_P]), ...
%     subs(4.0/sigma_g,[a,b,c,d, sigma_a, sigma_b, sigma_c, sigma_d], ...
%     [mu_R,mu_R,mu_R,mu_P, sigma_R, sigma_R, sigma_R, sigma_P]))


