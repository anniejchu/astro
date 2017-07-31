from visual import *
import numpy as np

scene.width = 900
scene.height = 600
scene.title = 'LENS EQUATION'

#SHOULD ONLY BE ABLE TO ENTER DL, SX,SZ, OZ,OX
ox = -80.0
oy = 0.0
oz = -40.0
opos = vector(ox, oy, oz)

lx = 0.0
ly = 0.0
lz = oz
lpos = vector(lx, ly, lz)


sx = 80.0
sy = 0.0
sz = 40.0
spos = vector(sx, sy, sz)

slmeetpos = vector(lx, ly, sz)

ds = sx-ox
dl = lx-ox
dls = sx-lx

ols = points(pos=[spos, lpos, opos, slmeetpos], size = 10, color = color.red)

#LINE FROM OBSERVER TO SOURCE
osline = paths.line(start=spos, end =opos, np = 20)
curve(pos=osline.pos)

#LINE FROM OBSERVER TO POINT WHERE SOURCE(parallel) AND LENS(perpendicular) MEET 
slmeetline = paths.line(start=opos, end=slmeetpos, np=20)
curve(pos=slmeetline.pos)

#LINE FROM OBSERVER TO SOURCE X AXIS
dsline = paths.line(start=opos, end=(sx, oy, oz), np = 20)
curve(pos=dsline.pos)

#LINE FROM X AXIS TO SOURCE
sperpline = paths.line(start=spos, end=(sx, oy, oz), np=20)
curve(pos=sperpline.pos)

#LINE FROM LENS TO SLMEET POINT (aka lens perp line)
slmeetlline = paths.line(start=lpos, end=slmeetpos, np = 50)
curve(pos=slmeetlline.pos)

#LINE FROM SLMEET POINT TO SOURCE (parallel)
slmeetsline = paths.line(start=slmeetpos, end=spos)
curve(pos=slmeetsline.pos)

#FINDING ANGLES 
ok = sz-oz #distance from x axis to source perpendicular
	#THETA
theta = atan2(ok,dl)
	#BETA
beta = atan2(ok, ds)
	#ALPHA
alpha = theta-beta

#beta opposite
oppbeta = tan(beta)*dl

#IMAGE STUFF
ok = sz-oz
iz = (tan(theta))*(ds)
ipos = vector(sx, sy, iz-sz)
#IMAGE LINE (FROM SLMEET TO IMAGE)
ps = points(pos=ipos, size = 10, color = color.cyan)
opsline = paths.line(start=slmeetpos, end=ipos, np = 20)
curve(pos=opsline.pos)

iperpline = paths.line(start=ipos, end=spos)
curve(pos=iperpline.pos)

#LABEL LINES
dsline = paths.line(start=(ox, oy, oz-20), end=(sx, oy, oz-20), np=20)
curve(pos=dsline.pos)
dlline = paths.line(start=(ox, oy, oz-10), end=(lx-1, oy, oz-10), np=20)
curve(pos=dlline.pos)
dlsline = paths.line(start=(lx+1, oy, oz-10), end=(sx, oy, oz-10), np=20)
curve(pos=dlsline.pos)


#LABELS
lenslabel = label(pos = lpos, text = 'lens', yoffset = -2, height = 10, color = color.white, line = False, box = False)
sourcelabel = label(pos = spos, text='source', yoffset = -2, height = 10, color = color.white, line = False, box = False)
observerlabel = label(pos = opos, text='observer', yoffset = -5, height = 10, color = color.white, line = False, box = False)
imagelabel = label(pos = ipos, text='image', yoffset = -2, height = 10, color = color.white, line = False, box = False)
slmeetlabel = label(pos= slmeetpos, text= 'sl meet', yoffset= -2, height = 10, color = color.white, line = False, box = False)

dlslabel = label(pos = ((sx+lx)*.5, oy, oz-15), text = 'DLS: '+ str(dls) +' pc', color = color.white, height = 8, line = False, box = False)
dllabel = label(pos = ((ox+lx)*.5, oy, oz-15), text = 'DL: '+ str(dl) +' pc', color = color.white,height = 8, line = False, box = False)
dslabel = label(pos = ((ox+sx)*.5, oy, oz-25), text = 'DS: '+ str(ds) +' pc', color = color.white, height = 8,line = False, box = False)

nlabel = label(pos = (sx, oy, (oz+sz)/2), text = 'n: '+str(sz-oz), color = color.cyan, height = 10, line = False, box = False)
zilabel = label(pos=slmeetpos, text = 'ξ: '+str((lz-oz)-oppbeta), yoffset = -slmeetpos.z, color = color.cyan, height = 10, line = False, box = False)


#show axes
arb = vector(150, 0, oz)
length = 50
sw = 0.5
	#X AXIS
xaxis = arrow(pos=(arb), axis = (1, 0, 0), length = length, color = color.white, shaftwidth = sw)
xalab = label(pos = (arb.x+length, arb.y, arb.z), text = 'x', color = color.white, line = False, box = False)
	#Y AXIS
yaxis = arrow(pos=(arb), axis = (0,1,0), length = length, color =color.white, shaftwidth = sw)
yalab = label(pos = (arb.x, arb.y+length, arb.z), text = 'y', color = color.white, line = False, box = False)
	#Z AXIS
zaxis = arrow(pos=(arb), axis = (0,0,1), length = length, color = color.white, shaftwidth = sw)
zalab = label(pos = (arb.x, arb.y, arb.z+length), text = 'z', color = color.white, line = False, box = False)


#show title
#title = label(pos = (-50, 0, 200), text = 'LENS EQUATION', height = 20, color = color.white)

#show angles
	#THETA
thetaarc = paths.arc(pos = opos, radius = 20, angle1 = 0, angle2 = -theta)
curve(pos=thetaarc.pos)
thetatext = label(pos = opos, text="θ", box = False, xoffset = -15, yoffset = 5, height=10)
	#BETA
betaarc = paths.arc(pos = opos, radius = 25, angle1 = 0, angle2 = -beta)
curve(pos=betaarc.pos)
betatext = label(pos = (opos.x+23, opos.y, opos.z+7), text = 'β', box = False, xoffset = -40, yoffset = 20, height = 10)
	#ALPHA
alphaarc = paths.arc(pos = opos, radius = 30, angle1 = -beta, angle2 = -theta)
curve(pos=alphaarc.pos)
alphatext = label(pos = (opos.x+20, opos.y, opos.z+17), text = 'α', box = False, xoffset = -45, yoffset = 20, height = 10)
	#DEFLECT ANGLE
deflectangle = paths.arc(pos = slmeetpos, radius = 20, angle1 = 0, angle2 = -theta)
curve(pos =deflectangle.pos)
deflecttext = label(pos = slmeetpos, text = "α'", box = False, xoffset = 20, yoffset = 10, height = 10) 

#print out data
ix = -120
iy = 0
iz = 100

init = vector(ix, iy, iz)

thetamas = 206265000*theta
betamas = 206265000*beta
alphamas = 206265000* alpha

#printout
print('theta (rad) : '+str(theta)+' (mas): '+str(thetamas))
print('beta (rad) : '+str(beta)+' (mas): '+str(betamas))
print('alpha (rad) : '+str(alpha)+' (mas): '+str(alphamas))