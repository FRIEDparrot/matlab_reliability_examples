clear, clc;
mu_ = [2.5, 2.5, 2.5, 2.5, 2.5];
sigma_d = [0.5, 0.5, 0.5, 0.5, 0.5];
sigma_ = diag(sigma_d.^2);
para_names = ["DS_FLAP_HEIGHT", "DS_LO_RIB_WIDTH1", "DS_LO_RIB_WIDTH2", ...
              "DS_LO_RIB_WIDTH3","DS_LA_RIB_WIDTH"];
g = @(x)2.5e8 - max_stress(x, para_names);   % 设定极限状态函数为最大应力不超过2.5e8Pa
initial_point = [2.0, 2.0, 2.0, 2.0, 2.0];

fx = @(x) joint_pdf(x, mu_, sigma_d);
params = MIAK_parainit(5);
params.num_MKIS_init = 60;
params.max_epoch1 = 20;

[dmodel, pref,msc]= MIAK_modelbuild(mu_, sigma_,g, initial_point,fx, params);

% 功能函数部分(需要传入各项的参数)
% 改良: 取消check_function, 直接在g(x)中处理传入的数据和进行数据检查
function yp = max_stress(xp, para_names)
    xp = round(xp, 4);    % retain for 4 digits of demical 
    m = size(xp,1);       
    yp = zeros (m, 1);   
    
    max_retry_time = 10;
    for i = 1:m
        tic;
        % 要求, 传入的xp中,每一项必须至少 >= 0.2, 否则直接认为失效点, 
        if min(xp(i,:), [], 2) < 0.2
            yp(i,:) = 6e8;        % 认为失效点的应力是6e8MPa
            fprintf("jumped invalid point %dth of total %d points, solution time: %f\n", i,m,  toc);
        else
            for retry = 1:max_retry_time
                f1 = fopen("journal.py", "r");
                cmd = string(fread(f1,'*char')');
                for j = 1:length(para_names)
                    cmd = strrep(cmd,para_names(j), string(num2str(xp(i,j))) );
                end
                % 将日志文件写入并且执行
                f2 = fopen("C:\Users\Parrot\Desktop\Recent\Graduate_Projects\Meta_IS_AK_Model\Meta_IS_AK_Test_files\scripts\journal.wbjn","w");
                fwrite(f2, cmd); 
                fclose all;
                system('"E:\Ansys2024R1\ANSYS Inc\v241\Framework\bin\Win64\RunWB2.exe" -B -R "C:\Users\Parrot\Desktop\Recent\Graduate_Projects\Meta_IS_AK_Model\Meta_IS_AK_Test_files\scripts\journal.wbjn"');
                
                Res = importdata("Stress_data.csv").data;
                xp_test = Res(size(Res, 1), 1:5);

                if find(xp_test ~= xp(i,:)) % 检查是否求解并导出结果成功
                    warning("Solution failed for design point %d\n");
                    if retry == max_retry_time
                        fprintf("The Solution terminated because of the ansys solution error\n");
                        error("Reach the Maximum Retrying %d, the solution failed\n", max_retry_time);
                    else
                        fprintf("         Retrying %d time of %d times.....\n",retry, max_retry_time);
                    end
                else
                    yp(i,:) = Res(size(Res, 1),6);
                    break;
                end
            end
            fprintf("finish solution point %d of total %d called, solution time: %f\n", i, m, toc);
        end
    end
end
