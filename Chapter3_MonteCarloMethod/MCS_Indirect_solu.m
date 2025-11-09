function [Pf,Pf_mu_X, Pf_D_X] = MCS_Indirect_solu(mu_,sigma_, g)
    %% %%%%%%%%%%%%%%%% 求解相关变量的部分 %%%%%%%%%%%%%%%%%%%%%%%% 
    num_MCS = 1e7;
    n = size(sigma_, 2);
    [A,D] = eig(sigma_); A = A';       % 分解为相互独立的随机变量, 然后转置矩阵;   D是y的sigma对应的矩阵(不是归一化相关系数矩阵)
    [~,~,v] = find(D);                 % 获取非零向量
    lambda_vector = v';                % 存储每一个特征值
    sigma_new_ = sqrt(lambda_vector);  % 转换后的方差
    
    % gy = @(y) max([g1((A'*y')'+ mu_), g2((A'*y')' + mu_)], [], 2);  % 也可以使用这个
    gy = @(y) g((A'*y')'+ mu_); %% ****** 对应公式3.4.2 ******* 
    yp = lhsnorm(zeros(1, n), D, num_MCS);  % 获取相应的MCS抽样点(y坐标下);
    
    Y = gy(yp);
    fail_points = find(Y < 0);
    
    Pf = size(fail_points, 1)/num_MCS;
    
    
    %% %%%%%%%%%%%%%% 求解正态空间中的可靠性灵敏度 (直接根据标准正态变量来求解 ) %%%%%%%%%%%%%%%%%
    sample_ = yp(fail_points,:);     %%%%%%  一定不要忘了整体取行
    % sample_res = Y(fail_points,:);   %%%%%%  对应的 y 部分
    
    % ------- 求解正态空间下的三个导数 -----------------
    Pf_D = zeros(n,n);  % 存储对于方差和相关系数的导数
    
    Pf_mu = sum((sample_ - 0)./(sigma_new_.^2))/num_MCS;   % 直接使用在正态空间中的部分求解均值
    Pf_sigma = sum(1./sigma_new_ .* ((sample_./sigma_new_).^2 -1))/num_MCS;  % Pf对simga导数
    for i = 1:n
        Pf_D(i,i) = Pf_sigma(i); 
    end
    % 按照3.4.4(3) 求解(其中前面的部分xTAx使用求和进行代替,后面加上均值);
    for i = 1:n
        for j = i+1:n
            [D_rho12, D_inv_rho12, D_det_rho12] = matrixdiff(D, i,j);  % 求解 Pf 对于 rho12的导数
            Pf_rho12 = -1/2 * ( sum(sum((sample_ * D_inv_rho12 .* sample_), 2))+ 1./det(D) * D_det_rho12)/num_MCS;
            Pf_D(i,j) = Pf_rho12; Pf_D(j,i) = Pf_rho12;
        end
    end
    % Pf_D 用于观察调试 
    %% %%%%%%%%%%%%% 转换到原变量空间下 (对应3.4.10)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [sigma_x, D_X] = nor_rhomat(sigma_);   % 求解x的方差和归一化相关系数矩阵D_X
    
    Pf_mu_X = sum(A .* Pf_mu, 2)';  % 注意是行相加的结果
    
    R_XY = zeros(n,n,n,n);   % RXY 是用于存储每个Y中的变量对X中的变量的导数 -> 以X对应元素下标为下标;
    temp_XY = zeros(n,n,n,n);       % 用于存储 sigma_ys, rho_ysr 对相应变量的导数
    
    for s = 1:n
        for r = 1:n
            a_s = A(s,:); a_r = A(r,:);
            if (s == r)
                % 对角线元素计算, sigma_Ys 对sigma_xi, rhoXiXj的偏导数
                temp = diag(sum(a_s' * a_s .* D_X .* sigma_x ./ (sigma_new_(s)), 2)); 
                for i = 1:n
                    for j = 1:n
                        if (i == j)
                            continue
                        end
                        temp(i,j) = 1/(2 * sigma_new_(s)) * a_s(i) * a_s(j) * sigma_x(i) * sigma_x(j);
                    end
                end
                temp_XY(:,:,s,s) = temp;   % sigma_Ys 对sigma_xi, rhoXiXj的偏导数
                clear a_s temp
            else % rho_Ysr 对 sigmaXi, rhoXiXj的导数
                % 求解对所有的 sigma_X 的导数
                temp = 2/(sigma_new_(s) * sigma_new_(r)) .* diag( sum(a_s' * a_r .* D_X .* sigma_x, 2));
                for i = 1:n
                    for j = 1:n
                        if (i == j)
                            continue
                        end
                        temp(i,j) = 1/(sigma_new_(s) * sigma_new_(r)) * a_s (i) * a_r(j) * sigma_x(i) * sigma_x(j);
                        % 系数不要写成2
                    end
                end
                temp_XY(:,:,s,r) = temp;
                clear a_s a_r temp
            end
        end
    end
    % temp_XY % 用于调试观察
    % 计算 R_XY 部分
    for i = 1:n
        for j = 1:n
            for u = 1:n
                for v = 1:n
                    R_XY(u,v,i,j) = temp_XY(i,j,u,v);
                end
            end
        end
    end
    clear temp_XY 
    
    Pf_D_X = zeros(n,n);
    for i = 1:n
        for j = 1:n
            Pf_D_X(i,j) = sum(sum(Pf_D .* R_XY(:,:,i,j)));
        end
    end
    %     Pf_mu_X
    %     Pf_D_X % 每个X对对应变量的导数
end

