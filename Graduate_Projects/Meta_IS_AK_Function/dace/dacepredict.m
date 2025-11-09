function [mu_g,sigma_g] = dacepredict(x, dmodel)
% 获取预测的均值和方差, 对于一个仍然适用;
    if size(x,1) == 1
		[mu_g, ~, sigma_g] = predictor(x,dmodel);
	else
		[mu_g, sigma_g] = predictor(x, dmodel);
    end
end