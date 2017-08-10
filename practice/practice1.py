#learning matplotlib.animation and vpython

import numpy as np 
from matplotlib import pyplot as plt 
from matplotlib import animation

#setting up the figure, the axis, and the plot element 
#we want to animate
fig = plt.figure() #creates figure window
ax = plt.axes(xlim=(0,2), ylim=(-2,2))
line, = ax.plot([], [], lw=2)

#initial function for plotting the background of every frame
def init():
	line.set_data([], [])
	return line,

#animation function
def animate(i): #only takes i as a single parameter ALWAYS i = frame number
	x = np.linspace(0,2,1000)
	y = np.cos(2*np.pi* (x-0.01*i))
	line.set_data(x,y)
	return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=50, interval=20, blit=True)
#fig is like void setup in processing
#animate changes the x and y every frame
#interval = delay in ms
#frames = how many frames it plays before restarting
#blit tells the animation to only redraw the pieces of the plot that changed

plt.show()
