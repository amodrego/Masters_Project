import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from os import listdir, remove
from sys import argv

import phisypy as phypy
import imageio  


FILE = argv[1]
SUBSTANCE = str(argv[2])

BASE_FILE = Path(FILE)
#BASE_FILE = Path('./output/')

# Takes the path given and runs the program with the files included in the folder
# BASE_FILE = Path(argv)
gif_name = 'output.gif' #BASE_FILE.name
counter = 0
ax = []
timestep = 120

substance_data = phypy.get_me_data(timestep, SUBSTANCE, BASE_FILE)
#colorbar = fig.add_axes([.92, .15, .012, .7])

sns.heatmap(substance_data, cmap="GnBu_r", linewidths=0.01, xticklabels=False, yticklabels=False)
filename = 'step_{}.png'.format(timestep)
ax.append(filename)
plt.savefig(BASE_FILE / filename)
plt.show()



        
