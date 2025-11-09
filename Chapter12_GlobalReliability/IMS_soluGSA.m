function GSA = IMS_soluGSA(mu_, sigma_, g, init_point, fx_pdf, hx_pdf, num_IMS)
%%%% 重要抽样法求解 %%%%%%%%%%%%%%%%%%%%%%%%% 
% fx_pdf : 原概率密度函数 
% hx_pdf : 重要抽样密度函数 
% init_point : 用于重要抽样的经验起始点, 
    n = size(mu_, 2);
    sigma_d = diag(sqrt(sigma_))';
    
    xp = slicesample(init_point,num_IMS,"pdf",hx_pdf);
    yp = g(xp);
    xp_fail = xp(yp < 0,:);
    
    Pf = mean((yp < 0) .* fx_pdf(xp) ./ hx_pdf(xp));  fprintf("Pf: %f\n", Pf);
    
    % 使用重要性采样方法, 将h(xi|F)投影到f(Xi|F)中, 获取类似于f(Xi|F)分布的失效样本点分布, 但不计算对应的功能函数
    fxF_samples = importance_sampling(hx_pdf, fx_pdf, xp_fail, size(xp_fail,1)); % f(Xi|F) 的 X 样本点
    
    %%%%%%%%%%%%%%%%%% 采用贝叶斯公式的思路, 求解在整个区域上的积分 %%%%%%%%%%%%%%%%%% 
    GSA = zeros(1,n);
    for i = 1:n
        [Pdf1, xi1] = ksdensity(fxF_samples(:,i)); % 使用核密度函数估计的策略, 直接生成f(X_i|F)分布
        
        fx_i = @(x) normpdf(x, mu_(i), sigma_d(i));           % fx(x)   的X概率密度函数
        fx_F = @(x) interp1(xi1,Pdf1,  x,"linear","extrap");  % fx(x|F) 的X概率密度函数
        GSA(i) = 1/2 * Pf * integral(@(x)abs(fx_i(x) - fx_F(x)), mu_(i) -  10 * sigma_d(i), mu_(i) + 10* sigma_d(i));
        fprintf("finish calculation of point %d/%d\n", i, n);
    end
end

% 重要性采样方法, 其中p(x)为抽样函数, q(x)为目标分布函数, 可以取得以 q(x) 的分布的样本
function q_samples = importance_sampling(p, q, p_samples, nums)
% 一般而言, 由于 p_samples 选取比 q_samples 多的值, 返回一个 nums 大小的向量
	w =  q(p_samples)./ p(p_samples);
	w = w/sum(w);
	idx = randsample(1:size(p_samples, 1), nums, true, w);
    q_samples = p_samples(idx,:);
end
