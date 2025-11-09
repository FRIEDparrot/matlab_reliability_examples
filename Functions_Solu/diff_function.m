function dg_dx = diff_function(g, x_i, delta)
    if nargin == 2
        delta = 1e-8; % 求导步进
    end
    % 求解一个函数g在x_i(行向量)处关于各个变量的导数, 注意输入x_i必须为行向量;
    
    m = size(x_i, 1); % 允许多行进行输入, 则返回多行的导数值
    n = size(x_i, 2);
    dg_dx = zeros(m, n);
    
    for i = 1:m
        cur_xi = x_i(i,:);
        for j = 1:n
            x_delta = zeros(1,n); x_delta(j) = delta;
            dg_dx(i, j) = (g(cur_xi + x_delta) - g(cur_xi))/delta;
        end
    end
end
