from visual import *

ball = sphere(pos = vector(-5,0,0), radius = 0.5, color = color.green)
wallR = box(pos = (6,0,0), size = (0.2, 12, 12), color = color.red)
wallL = box(pos = (-6, 0, 0), size = (0.2, 12, 12), color = color.red)
ball.velocity = vector(25, 5, 0)
vscale = 0.1
varr = arrow(pos=ball.pos, axis = vscale*ball.velocity, color = color.yellow)


deltat = 0.005
t=0


while t<3:
	if ball.pos.x > wallR.pos.x or ball.pos.x < wallL.pos.x:
		ball.velocity.x = -ball.velocity.x
		varr.axis.x = -varr.axis.x
	ball.pos = ball.pos + ball.velocity*deltat
	varr.pos = varr.pos + ball.velocity*deltat
#	arrow.axis = (vscale*ball.velocity) + arrow.pos
	t = t+deltat
	rate(100)


