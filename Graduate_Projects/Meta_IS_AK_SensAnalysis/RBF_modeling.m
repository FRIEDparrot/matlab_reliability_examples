%% %%%%%%%%%%%%%%% 基于线性响应面RBF方法, 实现RBF 模型的建立 %%%%%%%%%%%%%%%%
x = 0:0.01:4*pi;
y = sin(x);

load("E:\workpack\Matlab\MATLAB_reliability_engineering\Graduate_Projects\Meta_IS_AK_SensAnalysis\MIAK_results\IMS_result_final.mat")
net = newrb(x_train', y_train', 1e10);

sim(net, [0.2 0.2 0.2 0.2 0.2]')

sim(net,);