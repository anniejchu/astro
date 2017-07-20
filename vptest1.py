
from visual import *
r = vector(-3,4,0)

circle1 = sphere(pos = vector(1,1,1), radius = 0.25, color = color.cyan)
circle2 = sphere(pos = r, radius = 0.15, color = color.red)
pointer = arrow(pos = circle1.pos, axis=circle2.pos-circle1.pos, color = color.blue)


while circle2.x < 10:
	#sphere(pos=r, radius=0.5, color = color.cyan)
	rate(10)
	#circle2.pos = r
	circle2.x = circle2.x +1

