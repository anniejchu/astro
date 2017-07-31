from visual import *

scene.width = 600
scene.height = 600


for x in range(-400, 420, 50):
	for y in range(-400, 400, 50):
		sphere(pos=(x,y,0), radius = 5, color = color.white, opacity = 0.3)



