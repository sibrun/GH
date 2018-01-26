import matplotlib.pyplot as plt
import itertools
import logging


def save_matrix_plot(M, i_label, j_label, title, path):

    fig, ax = plt.subplots()
    ax.matshow(M, cmap=plt.cm.Blues)
    (i_max, j_max) = M.shape

    for (i, j) in itertools.product(range(i_max), range(j_max)):
        m = M[i][j]
        v = str(m) if m != -1 else ' '
        ax.text(j, i, v, va='center', ha='center')

    plt.title(title)
    plt.xlabel(j_label)
    plt.ylabel(i_label)
    plt.show()
    plt.savefig(path)
