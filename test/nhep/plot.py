import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


data_def = pd.read_csv("./default.csv")
data_15  = pd.read_csv("./p15.csv")
data_60  = pd.read_csv("./p60.csv")
data_120 = pd.read_csv("./p120.csv")
data_plog = pd.read_csv("./plog.csv")

labelsize = 13
legendsize = 13
ticksize = 12

plt.figure(figsize = (10,6))
plt.plot(data_def["timestep"]*1000, data_def["T"], label = "Default", color = "blue")
plt.plot(data_15["timestep"]*1000, data_15["T"], label = "Pref = 15 bar", color = "orange")
plt.plot(data_60["timestep"]*1000, data_60["T"], label = "Pref = 60 bar", color = "green")
plt.plot(data_120["timestep"]*1000, data_120["T"], label = "Pref = 120 bar", color = "red")
plt.plot(data_plog["timestep"]*1000, data_plog["T"], label = "Pressure-dependent (PLOG)", color = "black", marker = "o")
# for of plot a fake legend

plt.xlim(0.3,0.6)
plt.xlabel("Time (ms)", fontsize = labelsize)
plt.ylabel("T (K)", fontsize = labelsize)
plt.xticks(fontsize = ticksize)
plt.yticks(fontsize = ticksize)
plt.title("Comparison of OF and CTF Data", fontsize = labelsize)
plt.legend(frameon = False,loc = "best", fontsize = legendsize)
plt.tight_layout()
plt.savefig("../../results/p_comps.png", dpi = 300, bbox_inches = "tight")
