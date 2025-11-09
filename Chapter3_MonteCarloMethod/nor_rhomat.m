% 由sigma矩阵求解归一化相关系数矩阵
function [sigma_d, R] = nor_rhomat(D)
    sigma_d = sqrt(diag(D))';    % 求解每个sigma_;
    R = D./(sigma_d' * sigma_d);  % 求解归一化相关系数矩阵 rho
end
