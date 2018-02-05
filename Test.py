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

print([list(range(0,3)), list(range(3,6))])

L = SH.list_bipartite_g(2,3,4,3+2*5)

print(len(L))
for G in L:
    G.plot().show()








