# %%
import csv
import numpy as np
import matplotlib.pyplot as plt


def read_csv(csvpath, delimiter=','):
    f = open(csvpath, 'r')
    reader = csv.reader(f, delimiter=delimiter)
    index = 0
    header, datas = [], []
    for row in reader:
        if index == 0:
            temp = [item for number, item in enumerate(row)]
            header.append(temp)
        if index > 0:
            temp = [float(item) for number, item in enumerate(row)]
            datas.append(temp)
        index += 1
    header = np.array(header)
    datas = np.array(datas)
    f.close()
    return header, datas


def load_data_W_m2deq(filename):
    header, data = read_csv(filename)
    return data[:, 0], data[:, 4], data[0, 1]


# %% plot AR2_0deg data
plt.rcParams["text.usetex"] = True
plt.rcParams['font.family'] = 'serif'
plt.close('all')
fig, ax = plt.subplots(1, 1, figsize=(4, 3))


filename = "/media/gauthier/DataSSD/programs/gitLab/dry-granular/data/walls/AR2_0deg/AR2_0deg_3D.csv"

Ws, mu_eqs, H0 = load_data_W_m2deq(filename)
ax.plot((Ws*0.01/H0)**(-1), mu_eqs, 'p', c='cornflowerblue',
        label=r'3d B00 AR = 2  $\mu=0.44$', ms=4)
ax.errorbar((Ws*0.01/H0)**(-1), mu_eqs, yerr=0.02, c='cornflowerblue',
            fmt='none', markersize=8, markeredgewidth=1, capsize=2, )

# % plot AR0.5_15deg data
filename = "/media/gauthier/DataSSD/programs/gitLab/dry-granular/data/walls/AR0.5_15deg/AR0.5_15deg_3D.csv"

Ws, mu_eqs, H0 = load_data_W_m2deq(filename)

ax.plot((Ws*0.01/H0)**(-1), mu_eqs, 'o', c='darkred',
        label=r'3d B15 AR = 0.5  $\mu=0.44$', ms=4)
ax.errorbar((Ws*0.01/H0)**(-1), mu_eqs, yerr=0.02, c='darkred',
            fmt='none', markersize=8, markeredgewidth=1, capsize=2, )

ws_Ls=np.linspace(0.04,0.30,100)
muIonescu = 0.38 + 0.18 * 0.05 / ws_Ls
ax.plot((ws_Ls/H0)**(-1),muIonescu,'--',c='deeppink',label=' Ionescu et al. (2015) Correction',lw=1.4)
ax.set_position([0.15,0.15,0.8,0.8])

ax.set_ylabel(r'$\mu_{2D,eq}$')
ax.set_xlabel(r'$H_0/W$ ')
ax.legend()
plt.show()
