import GPy
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as ml
import matplotlib.patches as mpatches
import pandas as pd
import scipy.stats as stats

import time

X1 = np.array(pd.read_table("Rmatlab_files/generate_text/temp_to_X.txt", sep=","))
Y1 = np.array(pd.read_table("Rmatlab_files/generate_text/temp_to_r.txt", sep=",", header=None)[3])[:,None]
Y1
