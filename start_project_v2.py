
from visual import *
from visual.graph import *
import numpy as np
import pdb
#setup
scene = display(title = 'microlensing', width = 1000, height = 500, background = color.black)

scene.autoscale = False


#random variables
width = -300
height = 200
textsize = 10

#color variables 
red = (1,0,0)
green = (0,1,0)
blue = (0,0,1)
yellow = (1,1,0)
orange = (1,0.5, 0)
cyan = (0,1,1)
magenta = (1,0,1)
black = (0,0,0)
white = (1,1,1)
opacityplus = 0.4
opacityminus = 0.2
opacitycen = 0.3

#initial input variables 
imL = 10.0 #solar mass 
idL = 4000.0 #pc
idS = 8000.0 #pc
idLS = idS-idL
t0 = 57000.0
tr = 2000.0
muS =  np.array([8.0, 0.0])
muL =  np.array([0.00, 0.00])
beta = 1.8
t = np.arange(t0-tr, t0+tr)
x0S = 0.0
y0S = 0.0
x0L = -tr
y0L = 100.0
origin = vector(0,0, -idL)

# conversion
mL = imL * (1.99 * (10 ** 30))
dL = idL * (30856775714409184)
dS = idS * (30856775714409184)
dLS = dS-dL
G = 6.67 * 10**-11
c = 3.0 * 10**8
radtomas = 206265000
mtopc = 1/(3.0856776000000000 * 10**16)


#relative velocity vector
muRel = muS - muL
muRel1 = np.linalg.norm(muRel)

#calculating einstein radius 
inv_dist_diff = (1.0/dL)-(1.0/dS)
thetaE = (np.sqrt((4.0*G*mL/c**2) * inv_dist_diff)) # radians
thetaE1 = thetaE * radtomas #in mas
thetaE_hat = muRel/ thetaE1

#finding einstein radius distance from center
eradiusm = np.tan(thetaE)*dL #meters
eradiuspc = eradiusm*mtopc
eradiuspc_adjusted = eradiuspc*10**7


#einstein crossing time
tE = (thetaE1/muRel1) * 365 

#closest approach vector
u0_hat = np.zeros(2, dtype = float)
if beta > 0:
	u0_hat[0] = -np.abs(thetaE_hat[1])
	if np.sign(thetaE_hat).prod() > 0:
		u0_hat[1] = np.abs(thetaE_hat[0])
	else:
		u0_hat[1] = -np.abs(thetaE_hat[0])
else: 
	u0_hat[0] = np.abs(thetaE_hat[1])
	if np.sign(thetaE_hat).prod() > 0:
		u0_hat[1] = -np.abs(thetaE_hat[0])
	else:
		u0_hat[1] = np.abs(thetaE_hat[0])
u0_amp = beta / thetaE1
u0 = np.abs(u0_amp) * u0_hat

def getamp(u0=u0, thetaE_hat= thetaE_hat):
	tau = (t-t0)/ tE 
	u0 = u0.reshape(1, len(u0))
	thetaE_hat = thetaE_hat.reshape(1,len(thetaE_hat))
	tau = tau.reshape(len(tau), 1)

	#shape of u [n, 2]
	u = u0 + tau * thetaE_hat

	#shape of u amp
#	pdb.set_trace()

	u_amp = np.apply_along_axis(np.linalg.norm, 1, u)
	A = (u_amp**2 + 2)/(u_amp*np.sqrt(u_amp**2 +4))
	return A



#--------------DISPLAY---------------------------------------
#observers
OBPos = vector(0,0,1)
OB = sphere(pos = OBPos, radius = 1, color = white)

#SOURCE positioning
sPos = vector(x0S, y0S, -idL)+OBPos
SOURCE = sphere(pos = sPos, radius = 30, color=blue)

#LENS positioning 
lPos = vector(x0L, y0L, -idS)+OBPos
LENS = sphere(pos=lPos, radius = 50, color = yellow)
LENS.velocity = vector(1, 0, 0)
ER = ring(pos = lPos, radius = eradiuspc_adjusted, axis = (0,0,1), thickness=15, color= white)


thetaS = diff_angle(LENS.pos, SOURCE.pos)
lthetaplus = (thetaS+np.sqrt((thetaS**2)+(4*thetaE)))/2
lthetaminus = (thetaS-np.sqrt((thetaS**2)+(4*thetaE)))/2

thetaS1 = thetaS * radtomas
lthetaplus1 = lthetaplus* radtomas
lthetaminus1 = lthetaminus * radtomas

ldistplusm = np.tan(lthetaplus)*dLS
ldistpluspc = ldistplusm * mtopc
ldistpluspc_adjusted = ldistpluspc/(thetaE1*.8)

ldistminusm = np.tan(lthetaminus)*dLS
ldistminuspc = ldistminusm * mtopc
ldistminuspc_adjusted = ldistminuspc*10**5.5


cent_adjusted = (ldistpluspc_adjusted+ldistminuspc_adjusted)/2

#LETS CREATE THE LIGHT CURVES
plpos = vector(-ldistpluspc_adjusted,0,-idL)+OBPos
pluslight = sphere(pos=plpos, radius = 10, color = white, opacity = opacityplus)
plneg = vector(-ldistminuspc_adjusted, 0, -idL)+OBPos
minuslight = sphere(pos=plneg, radius = 10, color = white, opacity = opacityminus)
cenpos = vector(-cent_adjusted, 0, -idL)
cenlight = sphere(pos=cenpos, radius = 10, color = white, opacity = opacitycen)
#LABELS
lmasslabel = label(pos = origin, text = 'lens mass: '+ str(imL) +' solar masses', xoffset = -300, yoffset = 220, height = textsize, color = white, line = False)
ldistancelabel = label(pos = origin, text = 'distance to lens: '+str(idL)+' parsecs', xoffset = -280,yoffset= 195, height = textsize, color = white, line = False)
sdistancelabel = label(pos = origin, text = 'distance to source: '+str(idS)+ ' parsecs', xoffset = -266, yoffset = 170, height = textsize, color = white, line = False)
erlabel1 = label(pos = origin, text = 'einstein radius (angular): '+str(thetaE1)+ ' MAS', xoffset = -200, yoffset = 145, height = textsize, color =white, line = False)

llabel = label(pos = LENS.pos, text = 'L', height = textsize-2, color =white, line = False)
slabel = label(pos = SOURCE.pos, text = 'S', yoffset = 2, height = textsize, color =white, line = False)
erlabel = label(pos = ER.pos, text = 'ER', yoffset = llabel.yoffset+80, height = textsize-2, color = white, line = False)

amplabel = label(pos = origin, yoffset = -200, text = '', height = textsize, line=False)
tlabel = label(pos=origin, yoffset= -100, text = '', height = textsize, line = False)
sourceposlabel = label(pos = origin, yoffset = -200, xoffset = 300, text = '', height = textsize, line=False)
lensposlabel = label(pos = origin, yoffset = -175, xoffset = 300, text = '', height = textsize, line=False)



def moving(t = t, t0=t0, A = getamp()):
	lvel = LENS.velocity
	x = 0
	A1 = A**5
	#GRAPHS
	ampgraph = gdisplay(x=0, y = 300, width=500, height=300, title = 'AMP vs T', xtitle = 't', ytitle = 'amp', ymin = 1, ymax = A1[tr], xmin = t0-tr, xmax = t0+tr)
	ampcurve = gcurve(gdisplay = ampgraph, color = white)
	ampadjcurve = gcurve(gdisplay = ampgraph, color = white)

	lightcoordinatesx = gdisplay(x=600, y = 300, width=500, height=300, title = 'LIGHT CURVE(X) vs T', xtitle = 't', ytitle = 'x value',xmin = t0-tr, xmax = t0+tr)
	pluslightx = gcurve(gdisplay = lightcoordinatesx, color = cyan)
	minuslightx = gcurve(gdisplay = lightcoordinatesx, color = red)
	cenlightx = gcurve(gdisplay = lightcoordinatesx, color = white)

	lightcoordinatesy = gdisplay(x=600, y = 0, width=500, height=300, title = 'LIGHT CURVE(Y) vs T', xtitle = 't', ytitle = 'y value',xmin = t0-tr, xmax = t0+tr)
	pluslighty = gcurve(gdisplay = lightcoordinatesy, color = cyan)
	minuslighty = gcurve(gdisplay = lightcoordinatesy, color = red)
	cenlighty = gcurve(gdisplay = lightcoordinatesy, color = white)
	
	
	for time in t:
		LENS.pos = LENS.pos + lvel
		llabel.pos = llabel.pos + lvel
		ER.pos = ER.pos + lvel
		erlabel.pos = erlabel.pos + lvel
		xdiff = LENS.pos.x-SOURCE.pos.x
		if xdiff !=0:		
			rotateangle = np.arctan(200.0/xdiff)
		else:
			rotateangle = 0.0
		if time-t0 > -tr and time-t0 < 0.0:
			pluslight.pos.x = -np.cos(rotateangle)*ldistpluspc_adjusted
			pluslight.pos.y = -np.sin(rotateangle)*ldistpluspc_adjusted
			minuslight.pos.x = -np.cos(rotateangle)*ldistminuspc_adjusted
			minuslight.pos.y = -np.sin(rotateangle)*ldistminuspc_adjusted
			cenlight.pos.x = -np.cos(rotateangle)*cent_adjusted
			cenlight.pos.y = -np.sin(rotateangle)*cent_adjusted
		elif time-t0 > 0.0 and time-t0 < tr:
			pluslight.pos.x = np.cos(rotateangle)*ldistpluspc_adjusted
			pluslight.pos.y = np.sin(rotateangle)*ldistpluspc_adjusted
			minuslight.pos.x = np.cos(rotateangle)*ldistminuspc_adjusted
			minuslight.pos.y = np.sin(rotateangle)*ldistminuspc_adjusted
			cenlight.pos.x = np.cos(rotateangle)*cent_adjusted
			cenlight.pos.y = np.sin(rotateangle)*cent_adjusted


		pluslight.opacity = opacityplus*A1[x]
		minuslight.opacity = opacityminus*A1[x]
		cenlight.opacity = opacitycen*A1[x]

		amplabel.text = 'AMP: '+str(A1[x])+' (scaled e80)'+' REAL AMP: '+str(A[x])
		tlabel.text = 'Time: '+str(time-t0)	
		sourceposlabel.text = 'source x: '+str(SOURCE.pos.x)+' y: '+str(SOURCE.pos.y)
		lensposlabel.text = 'lens x: '+str(LENS.pos.x)+' y: '+str(LENS.pos.y)

		ampcurve.plot(pos=(time, A[x]))
		ampadjcurve.plot(pos=(time, A1[x]))
		pluslightx.plot(pos=(time, pluslight.pos.x))
		pluslighty.plot(pos=(time, pluslight.pos.y))
		minuslightx.plot(pos=(time, minuslight.pos.x))
		minuslighty.plot(pos=(time,minuslight.pos.y))
		cenlightx.plot(pos=(time, cenlight.pos.x))
		cenlighty.plot(pos=(time, cenlight.pos.y))


		x = x+1
		
		rate(500)
		

#graphamp()
moving()


