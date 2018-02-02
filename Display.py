import matplotlib.pyplot as plt
import numpy as np
import itertools


def save_2_indices_plot(value_dict, x_label, x_range, y_label, y_range, path):

    x_min = min(x_range)
    x_max = max(x_range)
    x_size = (x_max + 1 - x_min)*0.7
    y_min = min(y_range)
    y_max = max(y_range)
    y_size = (y_max + 1 - y_min)*0.7

    fig, ax = plt.subplots(figsize=(x_size, y_size))

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)

    for (x, y) in itertools.product(x_range, y_range):
        v = value_dict.get((x, y))
        if v is not None:
            ax.text(x, y, str(v), va='center', ha='center')

    x_ticks_grid = np.arange(x_min-0.5, x_max+1, 1)
    y_ticks_grid = np.arange(y_min-0.5, y_max+1, 1)

    x_ticks = np.arange(x_min, x_max+1, 1)
    y_ticks = np.arange(y_min, y_max+1, 1)

    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xticks(x_ticks_grid, minor=True)
    ax.set_yticks(y_ticks_grid, minor=True)

    ax.grid(which='minor')

    plt.tight_layout()

    plt.savefig(path)
