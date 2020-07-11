import sys
import json
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

from glob import glob

# ------------------------------------------------------------------------
# globals

colors = [
    '#79ADDC',
    '#FFC09F',
]

bins =  list(np.logspace(1.71, 3.95, 17))
bins += [1e10]

# ------------------------------------------------------------------------
# functions

def rm_prefix(s, p):
    if s.startswith(p):
        return s[len(p):]
    return s

def params(sim):
    sim  = sim.split('/')[0]
    p    = sim.split('_')[1:]
    p[0] = int(p[0]) # N
    p[1] = int(p[1]) # T
    # remainder are rates
    for i in range(2, len(p)):
        p[i] = float(p[i])
    return p

def coeffs(s):
    assert s[0] == '(' and s[-1] == ')'
    s = s[1:-1]
    p = s.split(',')
    return float(p[0]), float(p[1])

def adjust_spines(ax, spines, offset):
    """
    This is mostly from
    https://matplotlib.org/examples/pylab_examples/spine_placement_demo.html
    """
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', offset))  # outward by offset points
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])

# ------------------------------------------------------------------------
# Exploratory

# simple 1d histogram w/ length
def histogram_1d(d):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for i, (key, val) in enumerate(d.items()):
        mu, beta = coeffs(key)
        ax.hist(val["length"], bins=bins, color=colors[i], edgecolor='w', label=f"mu={int(mu):,}", density=True, alpha=.8)

    ax.set_xlabel("Length of block", fontsize='x-large')
    ax.set_ylabel('Fraction', fontsize='x-large')
    ax.set_xscale('log')
    ax.set_xlim([5e1, 9000])
    ax.tick_params(axis='both', which='major', labelsize='large')

    fig.legend(loc=(.75, .8), fontsize='large')
    adjust_spines(ax, ['left', 'bottom'], 5)

    fig.show()

def kde_2d(d):
    colors = ['Reds', 'Greens', 'Blues', 'Oranges', 'Purples']
    for i, (key, val) in enumerate(d.items()):
        mu, beta = coeffs(key)
        x, y = np.array(val['length']), np.array(val['depth'])
        ax = sns.kdeplot(data=x, data2=y, cmap=colors[i], shade=False, shade_lowest=False, alpha=.375, label=f"{int(mu):,}")

    ax.set_xlim([0, 10000])
    ax.set_xlabel("Block length (bp)", fontsize='x-large')
    ax.set_ylim([0, 1.25*N])
    ax.set_ylabel('Block depth (population count)', fontsize='x-large')

    ax.tick_params(axis='both', which='major', labelsize='large')
    leg = plt.legend(title='mu', title_fontsize='x-large', fontsize='large', markerscale=0.5)
    for lh in leg.legendHandles:
        lh.set_alpha(0.5)

# ------------------------------------------------------------------------
# Main point of entry

def as_array(*args):
    return [np.array(arg) for arg in args]

# TODO: better way to account/merge different runs...
def main(root):
    sims = glob(f"{root}/*/algo_stats.json")
    for sim in sims:
        N, T, mu, r_hgt, r_indel, r_xpose = params(rm_prefix(sim, root))
        # NOTE: this is a temporary hack
        if not (N == 100 and mu == 3e-4):
            continue

        with open(sim) as fd:
            d = json.load(fd)

        mus, mean_lens, sdev_lens, mean_deps, sdev_deps = [], [], [], [], []
        for i, (key, val) in enumerate(d.items()):
            mu, beta = coeffs(key)
            x, y = np.array(val['length']), np.array(val['depth'])
            mus.append(mu)
            mean_lens.append(np.mean(x))
            sdev_lens.append(np.std(x))
            mean_deps.append(np.mean(y))
            sdev_deps.append(np.std(y))
        mean_lens, sdev_lens, mean_deps, sdev_deps = as_array(mean_lens, sdev_lens, mean_deps, sdev_deps)

        idx = sorted(range(len(mus)), key= lambda i: mus[i])
        mean_lens = mean_lens[idx]
        sdev_lens = .5*sdev_lens[idx]
        mean_deps = mean_deps[idx]
        sdev_deps = .5*sdev_deps[idx]

        plt.errorbar(mean_lens, mean_deps, xerr=sdev_lens, yerr=sdev_deps)
        plt.xlabel("length")
        plt.ylabel("depth")
        plt.show()


if __name__ == "__main__":
    plt.style.use('seaborn-ticks')
    sys.exit(main(sys.argv[1]))
