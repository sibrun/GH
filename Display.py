import matplotlib.pyplot as plt
import numpy as np
import itertools

x_wieth = 0.7
y_width = 0.7
zero_symbol = '*'


def plot_2d_array(value_dict, x_label, x_range, y_label, y_range):

    x_min = min(x_range)
    x_max = max(x_range)
    x_size = (x_max + 1 - x_min) * x_wieth
    y_min = min(y_range)
    y_max = max(y_range)
    y_size = (y_max + 1 - y_min) * y_width

    fig, ax = plt.subplots(figsize=(x_size, y_size))

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)

    for (x, y) in itertools.product(x_range, y_range):
        v = value_dict.get((x, y))
        if v is not None:
            v = zero_symbol if v == 0 else v
            ax.text(x, y, str(v), va='center', ha='center')

    x_ticks_grid = np.arange(x_min - 0.5, x_max + 1, 1)
    y_ticks_grid = np.arange(y_min - 0.5, y_max + 1, 1)

    x_ticks = np.arange(x_min, x_max + 1, 1)
    y_ticks = np.arange(y_min, y_max + 1, 1)

    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xticks(x_ticks_grid, minor=True)
    ax.set_yticks(y_ticks_grid, minor=True)
    ax.grid(which='minor')

    plt.tight_layout()

    return plt


def plot_3d_array(value_dict, x_label, x_range, y_label, y_range, z_label, z_range, x_plots=2):

    x_min = min(x_range)
    x_max = max(x_range)
    x_size = (x_max + 1 - x_min) * x_wieth
    y_min = min(y_range)
    y_max = max(y_range)
    y_size = (y_max + 1 - y_min) * y_width
    z_min = min(z_range)
    z_max = max(z_range)
    z_size = z_max + 1 - z_min

    if z_size % x_plots:
        y_plots = int(round(float(z_size)/x_plots + 0.5))
    else:
        y_plots = z_size / x_plots

    fig, axarr = plt.subplots(y_plots, x_plots, figsize=(x_plots*x_size, y_plots*y_size))

    for z in z_range:
        ax = get_ax(axarr, x_plots, y_plots, z, z_min)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(str(z)+' '+z_label)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

    for (x, y, z) in itertools.product(x_range, y_range, z_range):
        v = value_dict.get((x, y, z))
        if v is not None:
            v = zero_symbol if v == 0 else v
            get_ax(axarr, x_plots, y_plots, z, z_min).text(x, y, str(v), va='center', ha='center')

    x_ticks_grid = np.arange(x_min - 0.5, x_max + 1, 1)
    y_ticks_grid = np.arange(y_min - 0.5, y_max + 1, 1)

    x_ticks = np.arange(x_min, x_max + 1, 1)
    y_ticks = np.arange(y_min, y_max + 1, 1)

    for z in z_range:
        ax = get_ax(axarr, x_plots, y_plots, z, z_min)
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        ax.set_xticks(x_ticks_grid, minor=True)
        ax.set_yticks(y_ticks_grid, minor=True)
        ax.grid(which='minor')

    plt.tight_layout()

    if z_size % x_plots:
        for z in range(z_max + 1, z_min + x_plots*y_plots):
            ax = get_ax(axarr, x_plots, y_plots, z, z_min)
            fig.delaxes(ax)

    return plt


def get_ax(axarr, x_plots, y_plots, z, z_min):
    dz = z - z_min
    if x_plots > 1 and y_plots > 1:
        i = dz / x_plots
        j = dz % x_plots
        ax = axarr[i, j]
    else:
        ax = axarr[dz]
    return ax
