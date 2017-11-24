#!/usr/bin/python3
# Takes a solution coordinate format and generates an image for
# each layer

import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas
import numpy as np

def print_usage():
  print("{} <solution.csv>".format(sys.argv[0]))

def print_cube(filename):
  data = pandas.read_csv("solution.csv")

  coords = data.as_matrix()
  maxval = np.max(data, axis=0)

  fig, ax = plt.subplots(1,maxval.z+1, figsize=(16,4));
  for z in range(0, maxval[2]+1):
    layer = np.zeros(shape=(maxval.x+1, maxval.y+1),dtype=int)
    layer_coords=data[(data.z == z)]
    for index, row in layer_coords.iterrows():
      layer[row.x, row.y] = row.v

    ax[z].imshow(layer,vmin=1,vmax=maxval.v, cmap="hsv")

    #text portion
    x, y = np.meshgrid(np.arange(0,maxval.x+1), np.arange(0,maxval.y+1))
    for x_val, y_val in zip(x.flatten(), y.flatten()):
      ax[z].text(x_val, y_val, chr(64+layer[y_val,x_val]), va='center', ha='center')

    ax[z].set_title("Layer {}".format(z))
  img_file="solution_layers.png";
  fig.savefig(img_file, dpi=100, bbox_inches="tight")
  print("Printed solution to {}".format(img_file))

if len(sys.argv) != 2:
  print_usage()
  exit(1)
print_cube(sys.argv[1])

