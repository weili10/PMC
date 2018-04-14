# author: 
# wei li      wli10@student.unimelb.edu.au
# side lu     sidel@student.unimelb.edu.au

# date: 21/10/2017   

# the program is used to plot the result of heat distribution

import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == '__main__':
    filename = sys.argv[1]
    a = np.loadtxt(filename)
    plt.imshow(a)
    plt.show()
