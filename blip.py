from visual import *
import numpy as np 

scene.width = 600
scene.height = 600
scene.autoscale = False
scene.range = 100

def stargrid(xdiff, ydiff): 
	for x in range(-40, 40, xdiff):
		for y in range(-40, 40, ydiff):
			sphere(pos=(x,y,0), radius = 5, color = color.white, opacity = 0.3)




def lensmoving():
	sx = random.randint(-50,50)
	sy = random.randint(-50,50)
	source = sphere(pos=(sx,sy,0), radius = 1.5, color = color.white, opacity = 0.1)

	lx = -90
	ly = -90
	lens = sphere(pos = (lx, ly, 0), radius = 1.5, color = color.cyan)
	lensv = vector(sx,sy,0)

	movex = (sx-lx)/50
	movey = (sy-ly)/50
	#tr = int(np.sqrt((sx-lx)**2+(sy-ly)**2))
	
	while lens.x < 100 and lens.y > -100 and lens.x > -100 and lens.y < 100:
		lens.x = movex + lens.x
		lens.y = movey + lens.y

		if lens.x-source.x < 0 and lens.x-source.x > -10:
			source.opacity = source.opacity * 1.5
		elif lens.x-source.x > 0 and lens.x-source.x < 10:
			source.opacity = source.opacity / 1.5
		#elif lens.x-source.x == 0:
		#	source.opacity = 1


		rate(10)


def dist():
	numpy.sqrt((sx-lx)**2+(sy-ly)**2)

	tr = np.absolute(sx-lx)
	t = np.arange(-tr, tr)

	timelab = label(pos = (-50, 50, 0), text = 1, line = False)


	for time in t:
		lens.x = movex + lens.x
		lens.y = movey + lens.y
		timelab.text = str(time)

		if time < 0 and time > -20:
			source.opacity = source.opacity + 0.05
		elif time > 0 and time < 20:
			source.opacity = source.opacity - 0.05


lensmoving()
#stargrid(100,120)