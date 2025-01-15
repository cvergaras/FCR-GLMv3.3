#%%
from glmpy import plots
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import os
#%%

lake = plots.LakePlotter("output/lake.csv")
fig, ax = plt.subplots(figsize=(10, 5))
lake.lake_level(ax=ax)


#%%
fig, ax = plt.subplots(figsize=(10, 5))
nc1 = plots.NCProfile("output/output.nc")
out = nc1.plot_var(ax=ax, var="PHQ_EP_Calcite")#, min_diff=0.01)
col_bar = fig.colorbar(out)
col_bar.set_label("Temperature (°C)")
plt.show()



# %%
nc1 = plots.NCProfile("output/output.nc")
vars = nc1.get_vars()
print(vars)

# %%
Data_GLM = nc.Dataset("output/output.nc")
ph = Data_GLM.variables["PHQ_EP_Calcite"][:, :, 0, 0]
print(ph[9])
# print(ph.min(), ph.max())


# %%
