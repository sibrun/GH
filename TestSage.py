import Display

if __name__ == '__main__':

    x_label = 'x'
    y_label = 'y'
    z_label = 'z'

    x_range = range(4,8)
    y_range = range(4,10)
    z_range = range(5,10)

    v_dict = {(5,4,5):1,(7,9,6):-1, (5,8,9):0}

    Display.plot_3d_array(v_dict, x_label, x_range, y_label, y_range, z_label, z_range, x_plots=3).savefig('./test.png')










