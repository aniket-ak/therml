import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import pandas as pd

import h5py
file1 = "/Users/aniket/Documents/MarlinSim/04_testing/scenarios/Calcite_Tiger/Solution/file_1.csv__solution.sol"
f = h5py.File(file1, "r")

time_duration = len(f['solution'].keys())
keys = list(f['solution'].keys())
keys_ = [float(i) for i in keys]
keys_.sort()
sol_shape = np.array(list(f['solution'][keys[0]])).T.shape

solution = np.empty((tuple([time_duration] + list(sol_shape))))

for i,k in enumerate(keys_):
    solution[i,:] = np.array(list(f['solution'][str(k)])).T

plt.matshow(solution[-1,:,:,0])
plt.colorbar()
plt.show()
