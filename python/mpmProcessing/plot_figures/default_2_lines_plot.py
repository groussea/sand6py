
#%%
import sys
import opyf  # from opyflow library some rendering function may be employed
sys.path.append(
    './..')
sys.path.append(
    './../..')
import d6py
import sys
sys.path.append(
    './../../imageProcessing')

from Tools_collapses import mask_collapses, mask_collapses2
import matplotlib.pyplot as plt
from d6py.Tools import *
import sys
import numpy as np
L = 20
H = 14

N = 2  # lignes
M = 3 # colonnes
Y_lim = [-1.5, H]
X_lim = [-L-0.2, 100]
w_axs = 0.23
h_axs = 0.45
w_s = 0.02
w_s2 = 0.005
cmapg = opyf.custom_cmap.make_cmap_customized(Palette='green')

colors = [(33./255, 66./255, 99./255, 0.05),
        (1, 1, 0.3, 0.9), (0.8, 0, 0, 0.9), (0, 0, 0, 0.9)]
position = [0, 0.1, 0.5, 1]
cmap = opyf.make_cmap(colors, position)

# cmap.set_over(color='g')
cmap.set_under(alpha=0)
plt.rcParams['font.family'] = 'serif'
plt.rc('text', usetex=True)
# %matplotlib qt5

#%%

def figure_2lines_tamplate():
    plt.close('all')
    plt.ion()

    fig, axs = plt.subplots(N, M, dpi=142, figsize=(4.88, 2.7))
    for i in range(N):
        for j in range(M):
            axs[i, j].plot([-L, -L], [-10, 2 * H], '-k', linewidth=2)
            axs[i, j].set_anchor('SW')
            axs[i, j].set_xlim([X_lim[0], 14.5])
            axs[i, j].set_ylim(Y_lim)
            axs[i, j].set_aspect('equal')
            axs[i, j].set_yticklabels([])
            axs[i, j].set_xticklabels([])

    for i in range(N): 
        axs[i,-1].set_xlim([X_lim[0], 70])


    for i in range(N):
        for j in range(M):
            axs[i, j].plot([-200, 200], [0, 0], 'k', linewidth=0.5)
            axs[i, j].plot([0, 0], [-200, 200], 'k', linewidth=0.5)
            axs[i, j].text(-0.5, 0, 'O', fontsize=10, horizontalalignment='right',
                        verticalalignment='bottom', bbox=dict(boxstyle="round", fc="w", alpha=0.2))
            draw_sc = 0
            if i == 0 and j == 0:
                x_sc, y_sc = -15, 0
                lx_sc, ly_sc = 5, 5
                draw_ax_sc(axs[i, j],x_sc, y_sc, lx_sc, ly_sc)

    axs[0,0].remove()
    axs[0,1].remove()
    axs[0,2].remove()
    axs[1,0].set_xlim([X_lim[0], 5])
    axs[1,1].set_xlim([X_lim[0], 40])
    e1=0.05 # ecart entre les deux axes du dessus
    l2 = 0.9 #longeur de lâxe du dessous
    i1 = (l2-e1)/2 #longueur de la figure du dessus
    pl = (1-l2)/2 #position left
    axs[1,0].set_position([pl,0.62,i1,0.28])
    axs[1,2].set_position([pl,0.2,l2,0.28])
    axs[1,1].set_position([pl+i1+e1,0.62,0.7,0.28])

    [x2, y2, X2, Y2] = axs[1, 2].get_position().bounds
    [x1, y1, X1, Y1] = axs[1, 1].get_position().bounds
    plt.pause(0.01)
    axs[1,1].set_position([x2+X2-X1,y1,X1,Y1])
    [x, y, X, Y] = axs[1, 0].get_position().bounds
    plt.figtext(x+X/2, y+Y+0.045, r'$t=0.2$ s',
                fontsize=9, horizontalalignment='center')
    [x, y, X, Y] = axs[1, 1].get_position().bounds
    plt.figtext(x+X/2, y+Y+0.045, r'$t=0.6$ s',
                fontsize=9, horizontalalignment='center')
    [x, y, X, Y] = axs[1, 2].get_position().bounds
    plt.figtext(x+X/2, y+Y+0.045, r'$t_f$',
                fontsize=9, horizontalalignment='center')
    plt.show()
    
    return fig, axs
# %%

def draw_gravity_2_lines(axs,slope,fontsize_g):
    from matplotlib.patches import Circle, Wedge, Polygon, Arc, Rectangle
    x_sc, y_sc = 45, 12
    lx_sc, ly_sc = 10, 5
    vecgy = y_sc
    vecgx = x_sc
    lg = 8
    dxT = 0.2*lg
    dyT = 0.1 * lg

    slopeinrad =slope * 3.15 / 180
    axs[1, 2].quiver(vecgx, vecgy, lg*np.sin(slopeinrad), -lg*np.cos(slopeinrad),
                    width=0.004, linewidth=1, angles='xy', scale_units='xy', scale=1)
    axs[1, 2].plot([vecgx, vecgx], [vecgy, vecgy-1.1*lg], 'k', linewidth=1)

    arc = Arc([vecgx, vecgy], lg*1.4, lg*1.4, theta1=270,
            theta2=270+np.rad2deg(slopeinrad), linewidth=0.5, color='k')
    axs[1, 2].add_patch(arc)
    axs[1, 2].text(vecgx+lg*slopeinrad+dxT, vecgy-lg *
                np.cos(slopeinrad), r'$\vec{g}$', fontsize=fontsize_g)
    axs[1, 2].text(vecgx+lg*slopeinrad+dxT, vecgy-lg*np.cos(slopeinrad)+5*dyT,
                r'$i$='+format(slope, '1.0f')+r'$^\circ $', fontsize=fontsize_g)