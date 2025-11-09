% @brief buncher 抽样函数
function xp = buncher_sample(x_i, sigma_d, varargin)
%buncher_sample buncher 抽样, 以点x_i为中心, 获取 2n+1 个buncher设计点
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
