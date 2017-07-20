from visual import *

ball = sphere(pos = vector(-5,0,0), radius = 0.5, color = color.green)
wallR = box(pos = (6,0,0), size = (0.2, 12, 12), color = color.red)
wallL = box(pos = (-6, 0, 0), size = (0.2, 12, 12), color = color.red)
wallUP = box(pos = (0,6,0), size = (12, 0.2, 12), color = color.cyan)
wallD = box(pos = (0,-6, 0), size = (12,0.2,12), color = color.cyan)
wallB = box(pos=(0,0,-6), size = (12, 12, 0.2), color = color.blue)
#wallF = box(pos = (0,0,6), size = (12,12, 0.2))

ball.velocity = vector(25, 5, 6)
ball.trail = curve(color = ball.color) #not actually a curve, ordered list of points

vscale = 0.1
varr = arrow(pos=ball.pos, axis = vscale*ball.velocity, color = color.yellow)


deltat = 0.005
t=0
scene.autoscale = False

while t<10:
	if ball.pos.x > wallR.pos.x or ball.pos.x < wallL.pos.x:
		ball.velocity.x = -ball.velocity.x
		varr.axis.x = -varr.axis.x
	if ball.pos.y < wallD.pos.y or ball.pos.y > wallUP.pos.y:	
		ball.velocity.y = -ball.velocity.y
		varr.axis.y = -varr.axis.y
	if ball.pos.z > 6  or ball.pos.z < wallB.pos.z:
		ball.velocity.z = -ball.velocity.z
		varr.axis.z = -varr.axis.z

	ball.pos = ball.pos + ball.velocity*deltat
	varr.pos = varr.pos + ball.velocity*deltat
	ball.trail.append(pos=ball.pos)
	t = t+deltat
	rate(100)


