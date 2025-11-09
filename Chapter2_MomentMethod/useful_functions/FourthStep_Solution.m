% 三点估计结合四阶矩方法获取的可靠性与可靠性灵敏度 
function [res_theta1, res_theta2, Failure_Properity] = FourthStep_Solution(mu_, sigma_, skewness_, kurtosis_, g)
    n = size(mu_, 2); res_theta1 = zeros(1,n); res_theta2 = zeros(1,n);
    l = zeros(n,3); p = zeros(n, 3);
    sigma_ = sigma_ .^ 0.5;
    % 其中l为三点部分, p为三点对应的权值, 每个变量对应三个点
    %%注意: sigma 是矩阵输入，因此必须使用 (i,i) 进行索引; 
    for i = 1:n      % 计算三个点以及对应的权值
        temp = sqrt(4 * kurtosis_(i) - 3 * skewness_(i)^2);     % 多次计算
        l(i, 1) = mu_(i) - sigma_(i,i)/2 * (temp - skewness_(i)); 
        l(i, 2) = mu_(i); 
        l(i, 3) = mu_(i) + sigma_(i,i)/2 * (temp + skewness_(i));
    
        p(i, 1) = 1/2 * (1 + skewness_(i)/temp )/(kurtosis_(i) - skewness_(i)^2);
        p(i, 2) = 1 - 1/(kurtosis_(i) - skewness_(i)^2);
        p(i, 3) = 1/2 * (1 - skewness_(i)/temp )/(kurtosis_(i) - skewness_(i)^2);
    end 
    %%%%%%%%%%%%%%%%%% 多变量情况求解的可靠性估计部分 %%%%%%%%%%%%%%%%%%%%%
    K = cycle_Matrix(n);  % 获取下标循环矩阵(每一个a对应的第二个下标)
    % 存储四个中心矩的计算值
    res_mu = 0; res_sigma = 0; res_skewness = 0; res_kurtosis = 0; 
    % K 是一个3^n x n 大小的矩阵, 计算每次中的P(i)和cur_point(i) 
    P = ones(3^n,1); Point_Matrix = zeros(3^n, n);   % P 是每一行给出的系数部分;
    
    for i = 1:3^n    % 行数
        index_ = K(i,:);  % 存储每一个Xi循环过程对应的下标值 i1, i2, ...
        for X_i = 1:n         % 变量个数的部分
            P(i) = P(i) * p(X_i, index_(X_i)); % 用来计算 E(Y^k)前面系数部分, 每一行存储一个数值, 避免重复计算;
            Point_Matrix(i,X_i) = l(X_i, index_(X_i));
        end
        cur_points = Point_Matrix(i,:);
        res_mu = res_mu + P(i) * g(cur_points);  % 仅用来计算均值变化 
    end 
    % 直接获取其中的P和cur_point部分, 避免重复计算
    for i = 1:3^n
        cur_points = Point_Matrix(i,:);
        temp = g(cur_points) - res_mu;
        res_sigma    = res_sigma  + P(i) * temp^2;
        res_skewness = res_skewness + P(i) * temp^3;
        res_kurtosis = res_kurtosis + P(i) * temp^4;
    end
    res_sigma = sqrt(res_sigma);  % 上面得到的数值是平方 
    res_skewness = res_skewness./(res_sigma.^3);
    res_kurtosis = res_kurtosis./(res_sigma.^4);
    
    beta_2M = res_mu/res_sigma; % 获取二阶可靠度指标和四阶可靠度指标;
    beta_4M = (3 * (res_kurtosis - 1) * beta_2M + res_skewness * (beta_2M^2 - 1))/sqrt((5 * res_skewness^2 - 9 * res_kurtosis  + 9) * (1 - res_kurtosis)); 
    Failure_Properity = normcdf(-beta_4M);   % 获取对应的失效概率;
    %%%%%%%%%%%%%%%%%%%     可靠性灵敏度计算部分    %%%%%%%%%%%%%%%%%%%%%%%
    L_norm = [-1/2 * sqrt(12), 0 , 1/2 * sqrt(12)];    % 正态分布的l取值点部分;
    % 计算对应的偏导数
    alp1g_thetak = zeros(1,2 * n); alp2g_thetak = zeros(1,2 * n); % 初始化部分
    for i = 1:3^n
        cur_points = Point_Matrix(i,:);
        for X_i = 1:n
            index_ = K(i,X_i);  % 存储每一个Xi循环过程对应的下标值 i1, i2, ...
            % 求解对于a3, a4的偏导数
            delta_1 = L_norm(index_) * P(i) * g(cur_points) ./ sigma_(X_i, X_i);
            delta_2 = (L_norm(index_)^2 - 1) * P(i) * g(cur_points) ./sigma_(X_i, X_i);
            delta_3 = L_norm(index_) * P(i) * (g(cur_points)- res_mu)^2 ./ sigma_(X_i, X_i);
            delta_4 = ((L_norm(index_)^2 - 1) * P(i) * (g(cur_points) - res_mu)^2)./ sigma_(X_i, X_i);
            alp1g_thetak(X_i) = alp1g_thetak(X_i) + delta_1;
            alp1g_thetak(X_i + n) = alp1g_thetak(X_i + n) + delta_2;
            alp2g_thetak(X_i) = alp2g_thetak(X_i) + delta_3;
            alp2g_thetak(X_i + n) = alp2g_thetak(X_i + n) + delta_4;
        end
    end
    % 初始化， 将alp3g_thetak和alp4g_thetak初始化成左右两项的和, 之后加上中间项
    alp2g_thetak = alp2g_thetak ./ (2 * res_sigma);
    alp3g_thetak = -3 * res_skewness ./ res_sigma * alp2g_thetak - 3 /res_sigma * alp1g_thetak;
    alp4g_thetak = -4 * res_kurtosis ./ res_sigma * alp2g_thetak - 4 * alp1g_thetak * res_skewness/res_sigma;
    
    % 求解中间项
    for i = 1:3^n
        cur_points = Point_Matrix(i,:);
        for X_i = 1:n
            index_ = K(i,X_i);  % 存储每一个Xi循环过程对应的下标值 i1, i2, ...
            alp3g_thetak(X_i) = alp3g_thetak(X_i) + L_norm(index_) * P(i) * (g(cur_points) - res_mu).^3 ./ res_sigma.^3 ./ sigma_(X_i, X_i);
            alp3g_thetak(X_i + n) = alp3g_thetak(X_i + n) + (L_norm(index_)^2 - 1) * P(i) * (g(cur_points) - res_mu).^3 ./ res_sigma.^3 ./ sigma_(X_i, X_i);
            alp4g_thetak(X_i) = alp4g_thetak(X_i) + L_norm(index_) * P(i) * (g(cur_points) - res_mu).^4 ./ res_sigma.^4 ./ sigma_(X_i, X_i);
            alp4g_thetak(X_i + n) = alp4g_thetak(X_i + n) + (L_norm(index_)^2 - 1) * P(i) * (g(cur_points) - res_mu).^4 ./ res_sigma.^4 ./ sigma_(X_i, X_i);
        end
    end
    Pf_beta4M = - exp(-beta_4M ^2./2)/(2 * pi).^0.5;
    beta4M_beta2M = (3 * res_kurtosis + 2 * res_skewness * beta_2M - 3)/sqrt((9 * res_kurtosis - 5 * res_skewness^2 -9) *(res_kurtosis - 1));
    beta2M_alp1g = 1/res_sigma;
    beta2M_alp2g = -res_mu /res_sigma^2;
    beta4M_alp3g = (beta_2M^2 - 1)./sqrt((9 * res_kurtosis - 5*res_skewness^2 - 9) * (res_kurtosis  -1)) + ...
        5 * (3 * (res_kurtosis - 1) * beta_2M + res_skewness * (beta_2M^2 -1)) * (res_kurtosis -1) * res_skewness/(2 * ((9 * res_kurtosis- 5 * res_skewness^2 - 9) * (res_kurtosis - 1))^1.5);
    beta4M_alp4g =    3 * beta_2M /sqrt((9 * res_kurtosis - 5 * res_skewness^2 -9) * (res_kurtosis - 1)) - ...
            (3 * (res_kurtosis - 1) * beta_2M + res_skewness * (beta_2M^2 -1)) * (18 * res_kurtosis - 18 - 5*res_skewness^2)/(2 * ((9 * res_kurtosis- 5 * res_skewness^2 - 9) * (res_kurtosis - 1))^1.5);
    
    SFp = Pf_beta4M.*(beta4M_beta2M .* (beta2M_alp1g .* alp1g_thetak + beta2M_alp2g .* alp2g_thetak) + beta4M_alp3g * alp3g_thetak + beta4M_alp4g * alp4g_thetak);
    
    for i = 1:n
        res_theta1(i) = SFp(i);
        res_theta2(i) = SFp(n+i);
    end
end

%%%%%%%%%%%%%%%    用于生成下标循环矩阵    %%%%%%%%%%%%%%%%%%%%%%%%%
function grid_ = cycle_Matrix(n)
    [K{1:n}] = ndgrid(1:3);  % 为每一个对应的变量建立一个 3x3 矩阵 (其中, i_k  = 1,2,3); 
    grid_ = reshape(cat(n+1, K{:}), 3^n, n); % 构造  * n 的矩阵, 每一格是grid经过reshape之后的值 
end