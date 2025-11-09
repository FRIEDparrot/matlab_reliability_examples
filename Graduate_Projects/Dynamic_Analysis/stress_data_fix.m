data = readmatrix("data_result.csv");


x =  data(:,2);
y = data(:,4);
z = data(:,5);


% scatter3(x,y,z)
% 将调节片真实模型的载荷投影到调节片简化模型范围上
% x_sim = [];
% y_sim = [];
% z_sim = [];

F = scatteredInterpolant(x,y,z);

x_new = linspace(min(x), max(x),300);
y_new =  linspace(min(y), max(y),90);

[xq, yq] = meshgrid(x_new, y_new);
zq = F(xq, yq);
for i = 1:size(zq,2)
    for j = 1:size(zq,1)
        if (x_new(i) >= 450 - abs(y_new(j)))
            zq(j,i) = NaN;
        end
    end
end

surf(xq, yq, zq);
shading interp
colormap default
colorbar
title("Pressure(MPa) distribution of the nozzle flap")

