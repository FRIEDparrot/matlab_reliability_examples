clear;
d = 4;
g=@(x) x(:,4)-x(:,2).*x(:,1)./(2*x(:,3)); %极限状态方程
x_mu = [460 20 19 392]; x_sigma = [7 2.4 0.8 31.4]; %均值及标准差

delta = ones(1,d);
x0=x_mu; %取初值为均值

x_design = AFOSM_solu(x_mu, diag(x_sigma.^2), g); % 设计点

q = qrandstream('sobol',d,'Skip',1e2);
N_IS = 1e4;

hx_pdf = @(x) prod(normpdf(x, x_design, x_sigma), 2);   % 提供重要抽样密度函数
x_IS = slicesample(x_design, N_IS ,"pdf", hx_pdf);      

fx = 1; hx = 1;
for i = 1:d
    fx = fx.*normpdf(x_IS(:,i),x_mu(i),x_sigma(i));
    hx = hx.*normpdf(x_IS(:,i),x_design(i),x_sigma(i));
end

y_IS = g(x_IS);

IF_IS = y_IS < 0;
Pf_IS = sum(IF_IS.*fx./hx)/N_IS; % 估计失效概率
x_f = x_IS(IF_IS,:); % 失效样本
% 马尔科夫链模拟
z_mar = ones(size(x_f));
z_mar(1,:) = x_f(1,:);
for k = 1:size(x_f,1)-1
    fx_k1 = 1; hz_k = 1; fz_k= 1; hx_k1 = 1;
    for i = 1:d
        fx_k1 = fx_k1.*normpdf(x_f(k+1,i),x_mu(i),x_sigma(i));
        hz_k = hz_k.*normpdf(z_mar(k,i),x_design(i),x_sigma(i));
        fz_k = fz_k.*normpdf(z_mar(k,i),x_mu(i),x_sigma(i));
        hx_k1 = hx_k1.*normpdf(x_f(k+1,i),x_design(i),x_sigma(i));
    end
    ratio = (fx_k1*hz_k)/(fz_k*hx_k1); % 计算比值
    rand_value = rand;
    
    if min(1,ratio) > rand_value % 确定下一个样本
        z_mar(k+1,:) = x_f(k+1,:);
    else
        z_mar(k+1,:) = z_mar(k,:);
    end
end

for i = 1:d
    xi_F = z_mar(:,i);
    lower_xi = x_mu(i) - 5*x_sigma(i);
    upper_xi = x_mu(i) + 5*x_sigma(i);
    xi_centers = linspace(lower_xi,upper_xi,30);
    pdf_fxi = normpdf(xi_centers,x_mu(i),x_sigma(i));
    ni_F = hist(xi_F,xi_centers);
    pdf_fxi_F = ni_F/length(xi_F)/(xi_centers(2)-xi_centers(1));
    delta(i) = sum(abs((pdf_fxi_F - pdf_fxi)*(xi_centers(2)-xi_centers(1))))*Pf_IS*0.5;
end

function f = mydiff1( fun,x,dim ) % 求偏导函数
    % x: 自变量(向量), x=[x1,x2,x3,...]
    % dim: 对第几个变量求偏导
    if dim<1,error('dim should >=1'),end;
    h=0.00001;
    n=length(x);
    if dim>n,error('dim should <=%d',n),end;
    I=zeros(1,n);
    I(dim)=1;
    f=(-fun(x+2*h*I)+8*fun(x+h*I)-8*fun(x-h*I)+fun(x-2*h*I))/(12*h);
end



% 使用马尔科夫链模拟,将h(xi)投影到f(Xi)中, 获取类似于f(Xi)分布的样本点, 但不计算对应的功能函数
% x_i = xp(1,:);                 % 初始化MCMC抽样的模拟点
% xpp = zeros(num_IMGRE, n);     % 实际符合f(X_i|F)分布的样本点
% 
% % 为了将 hx_pdf 投影到 fx_pdf 中 -> 下一个点的 fx ./ hx 和这个点进行比较
% for i = 1:num_IMGRE
%     r1 = (fx_pdf(xp(i,:))./ hx_pdf(xp(i,: ))) ./ (fx_pdf(x_i)  ./ hx_pdf(x_i));
%     r2 = rand();
%     % 使用 r = \frac{f1}/{h1} ./ {f2}/{f2} 的值
%     if (r2 < min(r1, 1))  % 当比值较小, 则进行转移
% 		x_i = xp(i,:);
%     end
%     xpp(i,:) = x_i;
% end
% 最终, 所获取到的即为 fx 分布的点
