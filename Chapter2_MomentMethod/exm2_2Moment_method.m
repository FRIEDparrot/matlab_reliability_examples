clear,clc
%%%%%% 
% 某内压圆筒型容器材料为15MnV, 随机变量取为内径D, 内压强P, 壁厚t与屈服强度sigma_s
% 这些随机变量独立且均服从正态分布, 参数在下面给出
D_mu = 460; D_sig = 7;
P_mu = 20;  P_sig = 2.4;
t_mu = 19;  t_sig = 0.8;
sigma_s_mu = 392; sigma_s_sig = 31.4 ; 


syms s p d t ss ps ds ts
g(s,p,d,t) = s - p * d/(2 * t);

% ------- derivates for function g------- 
g_s_ = subs(diff(g(s,p,d,t),s), [s,p,d,t],[sigma_s_mu, P_mu, D_mu, t_mu]);
g_p_ = subs(diff(g(s,p,d,t),p), [s,p,d,t],[sigma_s_mu, P_mu, D_mu, t_mu]);
g_d_ = subs(diff(g(s,p,d,t),d), [s,p,d,t],[sigma_s_mu, P_mu, D_mu, t_mu]);
g_t_ = subs(diff(g(s,p,d,t),t), [s,p,d,t],[sigma_s_mu, P_mu, D_mu, t_mu]);

mu_g = g(s,p,d,t);
sigma_g = sqrt(g_s_^2 * ss^2 + g_p_^2 * ps^2 + g_d_^2 * ds^2 + g_t_^2 * ts^2);

mu_g_var = vpa(subs(mu_g,[s,p,d,t], [sigma_s_mu,P_mu, D_mu,t_mu]));
sigma_g_var = vpa(subs(sigma_g,[ss,ps,ds,ts], [sigma_s_sig,P_sig, D_sig,t_sig]));

beta = mu_g/sigma_g;
beta_var = subs(beta,[s,p,d,t, ss,ps,ds,ts], [sigma_s_mu,P_mu, D_mu,t_mu, sigma_s_sig,P_sig, D_sig,t_sig]);

sprintf("mu_g: %f, sigma_g: %f",mu_g_var,sigma_g_var)
sprintf("beta: %f , Pf :%f", beta_var, normcdf(-beta_var))

beta_t_ = diff(beta,t);
beta_d_ = diff(beta, d);
beta_p_ = diff(beta, p);
beta_s_ = diff(beta, s);

Pf_beta_ = -1/sqrt(2 * pi) * exp(-mu_g^2/(2 * sigma_g^2));


sen_s = subs(beta_s_ * Pf_beta_, [s,p,d,t, ss,ps,ds,ts], [sigma_s_mu,P_mu, D_mu,t_mu, sigma_s_sig,P_sig, D_sig,t_sig] );
sen_p = subs(beta_p_ * Pf_beta_, [s,p,d,t, ss,ps,ds,ts], [sigma_s_mu,P_mu, D_mu,t_mu, sigma_s_sig,P_sig, D_sig,t_sig] );
sen_d = subs(beta_d_ * Pf_beta_, [s,p,d,t, ss,ps,ds,ts], [sigma_s_mu,P_mu, D_mu,t_mu, sigma_s_sig,P_sig, D_sig,t_sig] );
sen_t = subs(beta_t_ * Pf_beta_, [s,p,d,t, ss,ps,ds,ts], [sigma_s_mu,P_mu, D_mu,t_mu, sigma_s_sig,P_sig, D_sig,t_sig] );

sprintf("sense: %0.8f, %0.8f, %0.8f, %0.8f", sen_d, sen_p, sen_t, sen_s)

