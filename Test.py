import Display
from sage.all import *
import Shared as SH

if __name__ == '__main__':

    '''
    x_label = 'x'
    y_label = 'y'
    z_label = 'z'

    x_range = range(4,8)
    y_range = range(4,10)
    z_range = range(5,10)

    v_dict = {(5,4,5):1,(7,9,6):0, (5,8,9):1000}

    Display.plot_3d_array(v_dict, x_label, x_range, y_label, y_range, z_label, z_range,'./log/test.png', x_plots=3)
    '''

n_vertices = 8
n_loops = 7
n_edges = n_vertices + n_loops - 1

nauty_string = "-cbl -d1  %d %d:%d" % (n_vertices, n_edges, n_edges)
l = list(graphs.nauty_geng(nauty_string))

print(len(l))

for i in range(min(5,len(l))):
    g = l[i]
    g.plot().show()

