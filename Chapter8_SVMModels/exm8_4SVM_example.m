%% %%%%%%%%%%%%% Support Vector Machine %%%%%%%%%%%%% 
clearvars, clc;
%% %%%%%%%%%%%% 非线性单自由度无阻尼振动系统的SVM建模 %%%%%%% 
%mu_ = [1, 1, 0.1 , 0.5, 1, 1];
%sigma_d = [0.05, 0.1, 0.01, 0.05, 0.2, 0.2];
%omega_0 = @(x) sqrt((x(:,2) + x(:,3)) ./  x(:,1));
%g = @(x) 3 .* x(:,4) - abs( 2.* X(:,5) ./ (x(:,1).* omega_0(x).^2) .* sin(omega_0(x).^2 .* x(:,6)./2));

%% %%%%%%%%%%%% 两杆支撑系统的可靠性分析与SVM建模 %%%%%%%%%%% 
mu_ =    [200e6, 47.75e3, 100e-3, 100e-3,   0.03,  0.018, pi/3];
sigma_d = [20e6, 3.90e3,   5e-3,   5e-3, 1.5e-3, 0.9e-3, 1.15 * pi ./ 180];
g = @(x) x(:,1) - 2 .* x (:,2) .* sqrt(x(:,3).^2 + (x(:,4)./2).^2) ./(pi .* (x(:,5).^2 - x(:,6).^2)) .* (sin(x(:,7 ))./x(:,3) + 2 .* cos(x(:,7))./x(:,4));
sigma_ = diag(sigma_d.^2);

num_SVM = 1500;   % 支持向量机分类算法

Xtrain= lhsnorm(mu_, sigma_, num_SVM);
Ytrain = sign(g(Xtrain) + eps);
Xval = lhsnorm(mu_, sigma_, 500); 
Yval = sign(g(Xval));

% %%%%%%%%%%%%%%% create the SVC metamodel %%%%%%%%%%%%%%% 
uqlab 

%% %%%%%%%%%%%%% 支持向量分类模型预测方法 %%%%%%%%%%%%%%%%%%
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'SVC';               % use this for build the SVC metamodel
MetaOpts.ExpDesign.X = Xtrain;          
MetaOpts.ExpDesign.Y = Ytrain;

MetaOpts.ValidationSet.X = Xval;
MetaOpts.ValidationSet.Y = Yval;

MetaOpts.Kernel.Family = 'Gaussian';   % Select the Gaussian kernel family:
MetaOpts.EstimMethod   = 'CV';           % Cross-Validation 
MetaOpts.Optim.Method  = 'GS';           % Use the cross Enthalpy for the optimization;
MetaOpts.QPSolver      = 'SMO';          % SMO 序列最小优化算法
MetaOpts.Penalization  = 'linear';       % Use the span leave-one-out (LOO)
%%%%%%%%%%%%%% 建立分类代理模型 %%%%%%%%%%%%%%%%%%%% 
% --------- 生成新的样本, 使用模型进行预测----------
model = uq_createModel(MetaOpts);
X_new = lhsnorm(mu_, sigma_, 1e6);
[Y_class, Y_new] = uq_evalModel(model, X_new);
Pf = size(find(Y_new < 0), 1) ./ 1e6;
% uq_print(model);

fprintf("Pf for SVC: %f\n", Pf);
clear

%% %%%%%%%%%%%%% 支持向量回归模型预测方法 %%%%%%%%%%%%%%%%%%% 
MetaOpts.Type = 'Metamodel';              
MetaOpts.MetaType = 'SVR';                % use this for build the SVR metamodel
MetaOpts.ExpDesign.X = Xtrain;            
MetaOpts.ExpDesign.Y = g(Xtrain);

% ----- following are optional --------------------
% Select a linear penalization: 
% MetaOpts.Kernel.Family = 'Gaussian';     % Select the Polynomial kernel family-> note: ** select Gaussian may error!
% MetaOpts.EstimMethod   = 'CV';           % Also can use 'CV' as cross validation 
% MetaOpts.Optim.Method  = 'GA';           % Use the cross Enthalpy for the optimization;
% MetaOpts.QPSolver      = 'SMO';          % SMO 序列最小优化算法
% MetaOpts.Loss          = 'l1-eps';       % Loss-Function Type;
model = uq_createModel(MetaOpts);

% uq_print(model)    % show the basic infomation of the model builded
%%%%%%%%%%%%%% 建立回归代理模型 %%%%%%%%%%%%%%%%%%%% 
% --------- 生成新的样本, 使用模型进行---------- 
X_new = lhsnorm(mu_, sigma_, 1e6);
Y_new = uq_evalModel(model, X_new);
Pf = size(find(Y_new < 0), 1) ./ 1e6;
fprintf("Pf for SVR: %f\n", Pf);

Pf_mcs = MCS_solu(mu_, sigma_, g, 5e6);
fprintf("Pf for MCS: %f\n", Pf_mcs);