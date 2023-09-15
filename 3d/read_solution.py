import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import pandas as pd

import h5py
file1 = "/Users/aniket/Documents/MarlinSim/03_code/therml/3d/solution.jld"
f = h5py.File(file1, "r")

time_duration = len(f['solution'].keys())
sol_shape = np.array(list(f['solution']['1'])).T.shape

solution = np.empty((tuple([time_duration] + list(sol_shape))))
solution.shape
for i in range(1,time_duration):
    solution[i,:] = np.array(list(f['solution'][str(i)])).T

plt.matshow(solution[-1,:,:,1])