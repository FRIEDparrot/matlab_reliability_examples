clear, clc;
%% --------------- example 1 ----------------------------------------- 
mu_ = [0, 0];  sigma_d = [1, 1]; sigma_ = diag(sigma_d.^2);
c = 2;
g1 = @(x) c - 1 - x(:,2) + exp(- x(:, 1).^2 ./ 10) + (x(:,1)./5).^4;
g2 = @(x) c.^2 ./ 2 - x(:,1) .* x(:,2);
g = @(x) min([g1(x), g2(x)],[],2);
initial_point  = [-4, 4];             % 提供初始样本点(失效域内即可)
[dmodel, pref,msc]= MIAK_modelbuild(mu_, sigma_,g, initial_point);
fprintf("Pf = %f", msc.Pf)

%% --------------- example 2 -----------------------------------------
mu_ = [2,2,2,1]; 
sigma_d = [0.2, 0.2, 0.2, 0.25]; 
sigma_ = diag(sigma_d.^2);

g1 = @(x) 2 .* x(:,1) + 2 .* x(:,3) - 4.5 .* x(:,4);
g2 = @(x) 2 .* x(:,1) + x(:,2) + x(:,3) - 4.5 .* x(:,4);
g3 = @(x) x(:,1) + x(:,2) + 2 .* x(:,3) - 4.5 .* x(:,4);
g4 = @(x) x(:,1) + 2.* x(:,2) + x(:,3) - 4.5 .* x(:,4);
g = @(x) min([g1(x), g2(x), g3(x), g4(x)],[],2);
fx = @(x) joint_pdf(x, mu_, sigma_d);
params = MIAK_parainit(4);
params.max_epoch1 = 20;
params.nKMeans    = 15;

initial_point  = [1, 1, 1, 1];           % 提供初始样本点
[dmodel, pref, msc]= MIAK_modelbuild(mu_, sigma_,g, initial_point, fx, params);
fprintf("Pf = %f", msc.Pf)
