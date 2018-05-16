
__all__ = ['plot_array', 'plot_list', 'plot_2d_array', 'plot_3d_array']

import matplotlib.pyplot as plt
import numpy as np
import itertools
import logging
import pandas
import StoreLoad
import Parameters
import Shared


def plot_array(value_dict, ordered_param_range_dict, path, to_html=False, to_csv=False, x_plots=2, parameter_order=None):
    if len(ordered_param_range_dict) == 2:
        if parameter_order is not None:
            plot_2d_array(value_dict, ordered_param_range_dict, path, parameter_order=parameter_order)
        else:
            plot_2d_array(value_dict, ordered_param_range_dict, path)
    elif len(ordered_param_range_dict) == 3:
        if parameter_order is not None:
            plot_3d_array(value_dict, ordered_param_range_dict, path, x_plots=x_plots, parameter_order=parameter_order)
        else:
            plot_3d_array(value_dict, ordered_param_range_dict, path, x_plots=x_plots)
    if to_html or to_csv:
        plot_list(value_dict, ordered_param_range_dict, path, to_html=to_html, to_csv=to_csv)


def plot_list(value_dict, ordered_param_range_dict, path, to_html=True, to_csv=False):
    StoreLoad.generate_path(path)
    data_list = []
    for (key, value) in value_dict.items():
        if value is None:
           value = ' '
        elif value == '*':
            value = Parameters.zero_v_symbol
        data_list.append(list(key) + [value])
    data_list.sort()
    columns = ordered_param_range_dict.keys()+['dimension']
    data_frame = pandas.DataFrame(data=data_list, columns=columns)
    if to_html:
        html_path = path + '.html'
        data_frame.to_html(html_path)
    if to_csv:
        csv_path = path + '.csv'
        data_frame.to_csv(csv_path)


def plot_2d_array(value_dict, ordered_param_range_dict, path, parameter_order=(0, 1)):
    path += '.png'
    if parameter_order in {(0, 1), (1, 0)}:
        (x_idx, y_idx) = parameter_order
    else:
        raise ValueError('invalid parameter order')
    inverse_order = tuple(Shared.Perm(list(parameter_order)).inverse())

    (x_label, x_range) = ordered_param_range_dict.items()[x_idx]
    (y_label, y_range) = ordered_param_range_dict.items()[y_idx]
    if len(list(x_range)) == 0 or len(list(y_range)) == 0:
        logging.warn('empty parameter range: nothing to plot')
        return

    x_min = min(x_range)
    x_max = max(x_range)
    x_size = (x_max + 1 - x_min) * Parameters.x_width
    y_min = min(y_range)
    y_max = max(y_range)
    y_size = (y_max + 1 - y_min) * Parameters.y_width

    fig, ax = plt.subplots(figsize=(x_size, y_size))

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)

    for coordinates in itertools.product(x_range, y_range):
        old_coordinates = tuple(coordinates[i] for i in inverse_order)
        v = value_dict.get(old_coordinates)
        if v is not None:
            if v == '*':
                v = Parameters.zero_v_symbol
            elif v == 0:
                v = Parameters.zero_symbol
            (x, y) = coordinates
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

    StoreLoad.generate_path(path)
    plt.savefig(path)


def plot_3d_array(value_dict, ordered_param_range_dict, path, parameter_order=(0, 1, 2), x_plots=2):
    path += '.png'
    if parameter_order in {(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)}:
        (x_idx, y_idx, z_idx) = parameter_order
    else:
        raise ValueError('invalid parameter order')
    inverse_order = tuple(Shared.Perm(list(parameter_order)).inverse())

    (x_label, x_range) = ordered_param_range_dict.items()[x_idx]
    (y_label, y_range) = ordered_param_range_dict.items()[y_idx]
    (z_label, z_range) = ordered_param_range_dict.items()[z_idx]
    if len(list(x_range)) == 0 or len(list(y_range)) == 0 or len(list(z_range)) == 0:
        logging.warn('empty parameter range: nothing to plot')
        return

    x_min = min(x_range)
    x_max = max(x_range)
    x_size = (x_max + 1 - x_min) * Parameters.x_width
    y_min = min(y_range)
    y_max = max(y_range)
    y_size = (y_max + 1 - y_min) * Parameters.y_width
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

    for coordinates in itertools.product(x_range, y_range, z_range):
        old_coordinates = tuple(coordinates[i] for i in inverse_order)
        v = value_dict.get(old_coordinates)
        if v is not None:
            if v == '*':
                v = Parameters.zero_v_symbol
            elif v == 0:
                v = Parameters.zero_symbol
            (x, y, z) = coordinates
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

    StoreLoad.generate_path(path)
    plt.savefig(path)


def get_ax(axarr, x_plots, y_plots, z, z_min):
    dz = z - z_min
    if x_plots > 1 and y_plots > 1:
        i = dz / x_plots
        j = dz % x_plots
        ax = axarr[i, j]
    else:
        ax = axarr[dz]
    return ax
