clear, clc;
% 下面为失效函数的定义, 定义为并联失效模式的情况，即取系统的功能函数为max(g1, g2); 
% 有两个变量, 其系数已经在下面给出, 注意有相关系数0.4;
g1 = @(x) 2 - x(:,2) + exp(-0.1 .* x(:,1).^2) + (0.2.* x(:,1)).^4;
g2 = @(x) 4.5 - x(:,1) .* x(:,2);
g = @(x) max(g1(x), g2(x));

mu_ = [0.85, 0];
sigma_ = [(1/3)^2, 0.4 * 1/3; 0.4 * 1/3, 1^2];   % 此处sigma_为协方差矩阵;

%% *************** 求解相关变量的部分 ******************
num_MCS = 1e7;
n = size(sigma_, 2); 
[A,D] = eig(sigma_); A = A';    % 分解为相互独立的随机变量, 然后转置矩阵;
[~,~,v] = find(D);          % find(X) returns a vector containing the linear indices of each nonzero element in array X.
lambda_vector = v';             % 存储每一个特征值

Y = lhsnorm(zeros(1,n), D, num_MCS);  % 获取新的样本点(注意采用新的mu值)
X = (A' * Y' + mu_')';   % 注意, Y的每n组数据拼合为一个向量, 因此转置之后再相乘; -> 加上旧的mu值;

fail_points = find(g(X) < 0);
point_sample_ = X(fail_points,:);      % 失效点

Pf = size(fail_points, 1)/num_MCS        % 失效概率求解
Pf_Cov = sqrt((1 - Pf)/(num_MCS-1)/Pf)   % 失效概率方差
%% *************** 灵敏度求解部分 *************************
C1 = inv(sigma_);
% 均值灵敏度(3.4.4);
Pf_mu_ = sum((sigma_ \ (point_sample_ - mu_)'), 2)'./ num_MCS
% 求解协方差阵一般表达式:

%% ********** 求解协方差矩阵逆矩阵对于变量的偏导数 ********
delta = 0.0001; % 用于求解
sigma_d = sqrt(diag(sigma_)); % 求解每个sigma_;
Cov_Mat = sigma_./(sigma_d * sigma_d');  % 求解相关系数rho的矩阵

% 后面略去, 具体参考 (实际只是求解矩阵导数代入3.4.4即可)
% syms sigma1 sigma2 row12;
% Cx=[sigma1^2 row12*sigma1*sigma2; row12*sigma1*sigma2 sigma2^2];%协方差矩阵
% C2=inv(Cx); %协方差阵的逆矩阵
% dC2dsigma=[diff(C2,sigma1) diff(C2,sigma2)];%协方差阵关于标准差的偏导数
% D=det(Cx); %Cx的行列式
% dDdsigma=[diff(D,sigma1) diff(D,sigma2)];%协方差矩阵的行列式关于方差的导数
% dC2drow=[diff(C2,row12)]; %协方差矩阵关于相关系数的导数
% dDdrow=[diff(D,row12)]; %协方差矩阵的行列式对相关系数的导数
% %求关于方差的灵敏度
% S=subs(dC2dsigma,{sigma1,sigma2,row12},{1/3 1 0.4});
% T=subs(dDdsigma/D,{sigma1,sigma2,row12},{1/3 1 0.4});
% t1=sum(((x-miux)'*S(:,1:2))'.*(x-miux));
% t2=sum(((x-miux)'*S(:,3:4))'.*(x-miux));
% t=[(t1+T(1))' (t2+T(2))'];
% dpdsigma=-(0.5*I*t)/N
% %求解关于相关系数的灵敏度
% A=subs(dC2drow,{sigma1,sigma2,row12},{1/3 1 0.4});
% B=subs(dDdrow/D,{sigma1,sigma2,row12},{1/3 1 0.4});
% s1=sum(((x-miux)'*A)'.*(x-miux));
% s=(s1+B)';
% dpdrow=-(0.5*I*s)/N

