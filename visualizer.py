# Plot scatter plot from input.txt and make a polygon using points in output.txt
# Usage: python visualizer.py
# Output: scatter plot and polygon

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
import matplotlib.path as path

# Read input.txt
input_file = open("input.txt", "r")
# ignore first line
input_file.readline()
points = []
for line in input_file:
    points.append([float(x) for x in line.split()])
points = np.array(points)

# Read output.txt
output_file = open("output_merged.txt", "r")
polygon = []
for line in output_file:
    polygon.append([float(x) for x in line.split()])
polygon = np.array(polygon)

# Plot polygon
fig, ax = plt.subplots()
polygon = np.array(polygon)
polygon = np.vstack((polygon, polygon[0]))
codes = np.ones(len(polygon), int) * path.Path.LINETO
codes[0] = path.Path.MOVETO
path = path.Path(polygon, codes)
patch = patches.PathPatch(path, facecolor="orange", lw=2)
ax.add_patch(patch)

# Plot scatter plot
ax.scatter(points[:, 0], points[:, 1])
ax.plot(polygon[:, 0], polygon[:, 1])

#plot a scatter plot with output points
ax.scatter(polygon[:, 0], polygon[:, 1], c='r', s=50)
plt.show()
