import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot_band_structure(fig, ax, path="", prefix="graphene", valance_band_index=3, e_min = -5.0, e_max=5.0, k_special=[]):
    # load band structure data from file
    print('loading files from ' + path)
    data = np.loadtxt(path)
    k_vals = np.unique(data[:,0]) 
    bands = []
    num_k = len(k_vals)
    num_band = len(data[data[:,0] == k_vals[0]])
    # plot bands 
    fermi = max(data[valance_band_index*num_k:(valance_band_index + 1)*num_k, 1])
    for i in range(num_band):
        ax.plot(k_vals, data[i*num_k:(i + 1)*num_k, 1] - fermi, 'k', linewidth=0.5)
    ax.plot([min(k_vals), max(k_vals)], [0.0,0.0], color='gray', lw=0.5, alpha=0.6)
    # plot high symmetry lines
    for i in k_special: 
        x1 = [i, i]
        x2 = [e_min, e_max]
        ax.plot(x1, x2, 'k--', lw=0.5, alpha=0.5)
    # annotate plot 
    k_labels = [r'$\Gamma$', r'$K$', r'$M$', r'$\Gamma$']
    ax.set_xlim([min(k_vals), max(k_vals)])
    ax.set_ylim([e_min, e_max])
    plt.xticks(k_special, k_labels)
    plt.ylabel(r'Energy (eV)')
    # set title     
    title = prefix.capitalize()
    plt.title(title)
    fig.set_size_inches(2.2, 3.36)
    fig.savefig(prefix+"-bands.pdf", bbox_inches='tight', pad_inches=0.02)
    return fermi

def plot_dos_projected(fig, ax, path="", prefix="graphene", fermi=0.0, e_min = -5.0, e_max=5.0):
    # load band structure data from file
    print('loading files from ' + path)
    data = np.loadtxt(path)
    e_vals = data[:,0] - fermi
    dos_vals = data[:,1]*2 # a factor of two to account for two atoms in the unit cell
    # plot dos
    plt.fill_betweenx(e_vals, dos_vals, 0.0, color='k', lw=0.5, alpha=0.6)
    ax.plot([0, 3.0], [0.0,0.0], color='gray', lw=0.5, alpha=0.6)
    ax.set_xlim([0.0, 3.0])
    plt.xticks([0, 1, 2, 3], [0, 1, 2, 3])
    ax.set_ylim([e_min, e_max])
    plt.title(r'states/eV')
    plt.ylabel(r'Energy (eV)')
    plt.yticks([], [])
    fig.set_size_inches(1.16, 3.36)
    fig.savefig(prefix+"-dos.pdf", bbox_inches='tight', pad_inches=0.02)

fig, ax = plt.subplots(1, 1)
fermi = plot_band_structure(fig, ax, path="out/graphene.dat.gnu", prefix="graphene", valance_band_index=3, e_min = -20.0, e_max = 10.0, k_special=[0.00000, 0.6667, 1.0000, 1.5774])

# fig, ax = plt.subplots(1, 1)
# plot_dos_projected(fig, ax, path="out/graphene.pdos_tot", prefix="graphene", fermi=fermi, e_min = -20.0, e_max=10.0)
