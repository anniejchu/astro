#learning matplotlib.animation and vpython

import numpy as np 
from matplotlib import pyplot as plt 
from matplotlib import animation

#setting up the figure, the axis, and the plot element 
#we want to animate
fig = plt.figure() #creates figure window
ax = plt.axes(xlim=(0,2), ylim=(-2,2))
line, = ax.plot([], [], lw=2)

plt.show()