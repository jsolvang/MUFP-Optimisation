import numpy as np
import subprocess
import os
import matplotlib as plt

def plot_hydroD_results(results):
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    ax1.plot(x, y)
    ax1.set_title('Sharing x per column, y per row')
    ax2.scatter(x, y)
    ax3.scatter(x, 2 * y ** 2 - 1, color='r')
    ax4.plot(x, 2 * y ** 2 - 1, color='r')