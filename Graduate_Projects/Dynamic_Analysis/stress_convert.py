import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import griddata

stress_data = pd.read_csv("data.csv")

stress_data_res = pd.DataFrame({'points':[], 'x':[], 'y':[], 'z':[], 'p':[]})

index = 1
coefficient = 0.1

for i in range(stress_data.shape[0]):
    x = stress_data.iloc[i, 0]
    y = stress_data.iloc[i, 1]
    z = stress_data.iloc[i, 2]
    pres = coefficient * stress_data.iloc[i, 3];
    z_new = z * 150 / 465
    x_new = x * (450 - 150 * z_new/150)  / (960 - 495 * z/465);
    if (z != 0):
        line1 = np.array([index, x_new, 0, z_new, pres])
        line2 = np.array([index+1, x_new, 0, -z_new, pres])
        stress_data_res.loc[len(stress_data_res)] = line1   # 添加行
        stress_data_res.loc[len(stress_data_res)] = line2
        index = index + 2
    else:
        line1 = np.array([index, x_new, 0, z_new, pres])
        stress_data_res.loc[len(stress_data_res)] = line1
        index = index + 1

stress_data_res.to_csv("./data_result.csv",sep = ',', index=False)

print(stress_data_res)