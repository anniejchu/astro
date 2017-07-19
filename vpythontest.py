from visual import *

ball = sphere(pos = vector(-5,0,0), radius = 0.5, color = color.green)
wallR = box(pos = (6,0,0), size = (0.2, 12, 12), color = color.red)
ball.velocity = vector(25, 0, 0)
deltat = 0.005
t=0
while t<3:
	ball.pos = ball.pos + ball.velocity*deltat
	t = t+deltat


