clear, clc;
mu_ = [2.5, 2.5, 2.5, 2.5, 2.5];
sigma_d = [0.5, 0.5, 0.5, 0.5, 0.5];
sigma_ = diag(sigma_d.^2);
para_names = ["DS_FLAP_HEIGHT", "DS_LO_RIB_WIDTH1", "DS_LO_RIB_WIDTH2", ...
              "DS_LO_RIB_WIDTH3","DS_LA_RIB_WIDTH"];
g = @(x)2.5e8 - max_stress(x, para_names);   % 设定极限状态函数为最大应力不超过2.5e8Pa
% 
% [Pf, dmodel] = RS_Linear_solu(mu_, sigma_,g);
% save
