function [Pf, model] = RS_Linear_solu(mu_, sigma_, g)
% model: the predicted result is b1 * 1 + b2 * x(1) + b3 * x(2)...... 
%        use model.g for predict results 
    n =  size(mu_,2);
    sigma_d = diag(sqrt(sigma_))';
    weight_func = @(gx) min(abs(gx) + 1e-3)./(abs(gx) + 1e-3);   % WRSM 权重函数, g越接近0越大
    
    beta_pre = 0;              % 初始迭代beta值
    x0 = mu_;  y0 = g(mu_);    % 用于插值的初始样本点
    x_i = mu_;                 % 初始时, 设计点组以mu_ 为中心点;

    for epoch = 1:2000
	    % 第一次迭代 (buncher设计点选择方法), 选用新的2n个设计点(共有m = 2n+1个)
        xp = buncher_sample(x_i, sigma_d, 'f', 1);
        
        Y = g(xp); % 此处调用原始函数g求解所有点的 Y (m  * 1)矩阵;
        W = diag(weight_func(Y));                       % (m *  m);
        A = [ones(2 * n + 1,1), xp];                    % (m *n+1); 的拼接矩阵, 存储每一项系数
        b = ((A' * W * A) \ (A' * W * Y))';             % 获取待定系数向量 (n+1 * 1)
        % 对应的响应面函数为b0 + b1 x1 + b2 x2 + ... + bn xn,
        
        % 将预测响应面函数代替原响应函数, 代入AFOSM中
        g_new = @(x) sum(b.* [ones(size(x,1), 1), x], 2);     % 返回响应面模型系数
        % 用 AFOSM 求解设计点的可靠度指标 beta^k 以及设计点 -> 实际上最后还是AFOSM求解, 对于高次的往往不适用
        [x_i, beta_res, Pf] = AFOSM_solu(mu_, sigma_, g_new); 
        
        model.b = b;
        model.g = @(x)sum(model.b.* [ones(size(x,1), 1), x], 2);
        % 分步保存
        % if mod(epoch, 5) == 0
            % save(join(["FunctionSolu\backups\RSL_data_epoch", num2str(epoch), ".mat"],''),"Pf","model","beta_res", "x_i","Y","b");
        % end
        fprintf("epoch: %d, Pf: %f, beta : %f, beta_variance: %f(target: 0.001)\n", epoch, Pf, beta_res, abs(beta_res - beta_pre));
        if abs(beta_res - beta_pre) < 0.001
            break;
        else  % 进行插值,  更新设计点
            x_i = mu_ +  (x_i - x0) * y0 / (y0 - g_new(x_i)); 
            beta_pre = beta_res;
        end
    end
end

% @brief buncher 抽样函数
% buncher_sample buncher 抽样, 以点x_i为中心, 获取 2n+1 个buncher设计点
function xp = buncher_sample(x_i, sigma_d, varargin)
%f 插值系数, 默认为1
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);   % 检查输入的部分是不是正整数
    addRequired (p,'x_i');      % 占据第一个位置
    addRequired (p,'sigma_d');  % 占据第二个位置 
    addOptional(p, 'f', 1, validScalarPosNum);  % 添加可选参数(插值系数f)
    parse(p, x_i, sigma_d ,varargin{:});          % 匹配输入参数, 同时进行检查

    f = p.Results.f;            % 获取对应的参数

    n = size(x_i,2);
    xp = zeros(2 * n + 1, n);   % 设计点组
    xp(1,:) = x_i;
    for i = 1:length(x_i)
        f_delta = zeros(1,n); f_delta(i) = sigma_d(i) * f;
        xp(2 * i,:) = x_i + f_delta;
        xp(2*i+1,:) = x_i - f_delta;
    end
end
