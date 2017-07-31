from visual import *
import numpy as np 

scene.width = 600
scene.height = 600
scene.autoscale = False
scene.range = 100
scene.title = 'BLIP'

def randostars():
	rx = []
	ry = []
	opacity = []
	radius=  []
	for x in range(0, 7):
		rx.append(random.randint(-90,90))
		ry.append(random.randint(-90,90))
		opacity.append(random.uniform(0.2, 0.8))
		radius.append(random.uniform(1,2))
		sphere(pos=(rx[x],ry[x],0), radius = radius[x], color = color.white, opacity = opacity[x])



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
	
	while lens.x < 100 and lens.y > -100 and lens.x > -100 and lens.y < 100:
		lens.x = movex + lens.x
		lens.y = movey + lens.y

		if lens.x-source.x < 0 and lens.x-source.x > -10:
			source.opacity = source.opacity * 1.5
		elif lens.x-source.x > 0 and lens.x-source.x < 10:
			source.opacity = source.opacity / 1.5


		rate(10)

randostars()
lensmoving()
