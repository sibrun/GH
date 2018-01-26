import numpy as np
import Display
import matplotlib.pyplot as plt
import os

reload(Display)

min_val, max_val = 0, 15

M = np.random.randint(0, 10, size=(max_val, max_val-1))

path = os.path.join('test','test.png')

Display.save_matrix_plot(M,'x','y','titel',path)

