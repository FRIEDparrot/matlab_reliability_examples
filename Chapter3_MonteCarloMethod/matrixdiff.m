%%%%%% 根据协方差矩阵获取矩阵和det对某个\rho或者sigma, 或者变量导数的代码 %%%%%%%
function [res_diff, res_invdiff, res_detdiff] = matrixdiff(D_X, row, col)
delta_x = 1e-9;  % 变量增量

if (length(row) ~= length(col))
    error("the row array and col array must be match");
end
n  = size(D_X, 1);
m = length(row);
sigma_d = sqrt(diag(D_X));            % 求解每个sigma_;
Cov_Mat = D_X./(sigma_d * sigma_d');  % 求解归一化相关系数矩阵 rho 
C_x_inv = inv(D_X); C_x_det = det(D_X);

res_diff = zeros(n,n,m); res_invdiff = zeros(n,n,m); res_detdiff = zeros(m);

for i = 1:m
    if row(i)~=col(i) % 求解相关系数的导数(对称增加);
        a = row(i);  b = col(i);
        delta_mat = zeros (n,n); delta_mat(a,b) = delta_x; delta_mat(b,a) = delta_x;
        Cov_Mat_new = Cov_Mat + delta_mat; % 增量
        sigma_d_new = sigma_d;
    else % 求解方差的导数
        a = row(i);
        delta_mat = zeros(1,n); delta_mat(i) = delta_x;
        sigma_d_new = sigma_d + delta_mat;
        Cov_Mat_new = Cov_Mat;
    end
    C_x_New = sigma_d_new * sigma_d_new' .* Cov_Mat_new;  % 注意此处使用.*
    C_x_inv_New = inv(C_x_New);
    C_x_det_New = det(C_x_New);
    
    res_diff(:,:,i) = (C_x_New - D_X)./ delta_x;
    res_invdiff(:,:,i) = (C_x_inv_New - C_x_inv) ./ delta_x;
    res_detdiff(i) = (C_x_det_New - C_x_det)/delta_x;
end
