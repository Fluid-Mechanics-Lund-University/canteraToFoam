import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data_of  = pd.read_csv("./of.csv")
data_ctf = pd.read_csv("./ctf.csv")
labelsize = 13
legendsize = 13
ticksize = 12

plt.figure(figsize = (10,6))
for i in range(len(data_of.columns)):
    if data_of.columns[i] == "timestep":
        continue
    max_of = data_of.iloc[:, i].max()
    max_ctf = data_ctf.iloc[:, i].max()
    of_plot  = data_of.iloc[:, i] / max_of
    ctf_plot = data_ctf.iloc[:, i] / max_ctf
    max_value_difference = np.abs(max_of - max_ctf)/(np.abs(max_of) + 1e-15)
    plt.plot(data_of["timestep"], of_plot)
    labelString = "CTF:" + data_of.columns[i].replace("_", " ")
    labelString += r" ($Diff_{max}$: " + f"{max_value_difference*100:.2f}%)" 
    plt.plot(data_ctf["timestep"], ctf_plot, linestyle='--', label = labelString)
# for of plot a fake legend
plt.plot([],[] , label = "OF")
plt.xlabel("Time (s)", fontsize = labelsize)
plt.ylabel("Normalized Value", fontsize = labelsize)
plt.xticks(fontsize = ticksize)
plt.yticks(fontsize = ticksize)
plt.title("Comparison of OF and CTF Data", fontsize = labelsize)
plt.legend(frameon = False,loc = "center left", fontsize = legendsize)
plt.tight_layout()
plt.savefig("../../results/ctf_of_comps.png", dpi = 300, bbox_inches = "tight")