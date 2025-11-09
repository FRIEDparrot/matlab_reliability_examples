% @bug: Pf 不能使用这种方法直接求解, 需要预测求解 
function [g_new, beta_res, Pf] = RSM_solu(mu_, sigma_ , g)
    % RSM_SOLU 线性加权响应面方法进行的可靠性分析
    n = size(mu_, 2);
    sigma_d = sqrt(diag(sigma_));
    %% %%%%%%%%%%%%%%%%%%% 利用计算权重 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    weight_func = @(gx) min(abs(gx) + 1e-3)./(abs(gx) + 1e-3);   % WRSM权重函数, g越接近0越大
    % weight_func = @(gx) ones(length(x), 1);  % TWRSM权重函数
    
    beta_pre = 0;              % 初始迭代beta值
    x0 = mu_;  y0 = g(mu_);    % 用于插值的初始样本点
    x_i = mu_;                 % 初始时, 设计点组以mu_ 为中心点;
    for epoch = 1:1000
        % 每一次迭代 (buncher设计点选择方法), 选用新的2n个设计点 
        xp = buncher_sample(x_i, sigma_d,'f',1);
        m = length(xp);  % 构造以下矩阵
        Y = g(xp);  % 此处调用原始函数g求解所有点的Y (m  * 1);
        W = diag(weight_func(Y));               % (m *  m);
        A = [ones(m,1), xp];                    % (m *n+1);
    
        b = ((A' * W * A) \ (A' * W * Y))';  % 获取待定系数向量 (n+1 * 1)
        % 对应的响应面函数为b0 + b1 x1 + b2 x2 + ... + bn xn,
        % 只需将预测响应面函数代入AFOSM中代替原响应函数即可
        
        % ------------- 构造预测函数 g(x) ---------------- 
        g_new = @(x) b(1);  
        for i = 1:n
            g_new =@(x) g_new(x) + b(i+1) .* x(:,i);
        end

        % 用 AFOSM 求解设计点的可靠度指标 beta^k 以及设计点
        [x_i, beta_res, Pf] = AFOSM_solu(mu_, sigma_, g_new); 
        
        if abs(beta_res - beta_pre)/beta_pre < 0.001
            break;
        else  % 进行插值,  更新设计点
            x_i = mu_ +  (x_i - x0) * y0 / (y0 - g(x_i));
            sprintf("epoch: %d, Pf: %f, beta : %f", epoch, Pf, beta_res)
            beta_pre = beta_res;
        end
    end
    
end