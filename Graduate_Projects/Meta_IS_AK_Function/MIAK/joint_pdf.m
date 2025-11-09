function f_x = joint_pdf(x, mu_,sigma_d)
% brief: 求解多变量的联合概率密度函数(normpdf)
    f_x = prod(1./sqrt(2 .* pi .* sigma_d.^2).* exp(-0.5 .* ((x - mu_)./(sigma_d)).^2), 2);
end
