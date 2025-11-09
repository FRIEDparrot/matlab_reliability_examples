%% %%%%%%%%%%%%%%%%%% 使用RBF方法建立神经网络,获取在相同验证集下的 %%%%%%%%%%%%%%%%%%% 

data = load("whitenoise_stress_data.csv");
t = data(:,1);
p = data(:,2);

plot(t,p)