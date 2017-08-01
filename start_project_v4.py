
from visual import *
from visual.graph import *
import numpy as np
import pdb
#setup
scene = display(title = 'microlensing', width = 750, height = 350, background = color.black)

scene.autoscale = False


#random variables
textsize = 10
xoff = -200


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



# conversions and constants

G = 6.67 * 10**-11
c = 3.0 * 10**8
radtomas = 206265000
mtopc = 1/(3.0856776000000000 * 10**16)
smtokg = 1.99 * (10 ** 30)
pctom = 30856775714409184

class PSPL(object):
	def __init__(self, imL, idL, idS, t0, tr, muS, muL, beta, x0S, y0S, y0L):
		self.imL = imL
		self.idL = idL
		self.idS = idS
		self.t0 = t0
		self.tr = tr
		self.muS = muS
		self.muL = muL
		self.beta = beta
		self.x0S = x0S
		self.y0S = x0S
		self.y0L = y0L

		#basic
		self.t = np.arange(self.t0-self.tr, self.t0+self.tr)
		self.mL = self.imL * smtokg
		self.dL = self.idL * pctom
		self.dS = self.idS * pctom
		self.dLS = self.dS-self.dL
		self.x0L = -self.tr
		#relative velocity vector
		self.muRel = self.muS - self.muL
		self.muRel1 = np.linalg.norm(self.muRel)

		#calculating einstein radius 
		inv_dist_diff = (1.0/self.dL)-(1.0/self.dS)
		self.thetaE = (np.sqrt((4.0*G*self.mL/c**2) * inv_dist_diff)) # radians
		self.thetaE1 = self.thetaE * radtomas #in mas
		self.thetaE_hat = self.muRel/ self.muRel1

		#finding einstein radius distance from center
		eradiusm = np.tan(self.thetaE)*self.dL #meters
		eradiuspc = eradiusm*mtopc #pc
		self.eradiuspc_adjusted = eradiuspc*10**7


		#einstein crossing time
		self.tE = (self.thetaE1/self.muRel1) * 365

		#closest approach vector
		self.u0_hat = np.zeros(2, dtype = float)
		if beta > 0:
			self.u0_hat[0] = -np.abs(self.thetaE_hat[1])
			if np.sign(self.thetaE_hat).prod() > 0:
				self.u0_hat[1] = np.abs(self.thetaE_hat[0])
			else:
				self.u0_hat[1] = -np.abs(self.thetaE_hat[0])
		else: 
			self.u0_hat[0] = np.abs(self.thetaE_hat[1])
			if np.sign(self.thetaE_hat).prod() > 0:
				self.u0_hat[1] = -np.abs(self.thetaE_hat[0])
			else:
				self.u0_hat[1] = np.abs(self.thetaE_hat[0])
		self.u0_amp = self.beta / self.thetaE1
		self.u0 = np.abs(self.u0_amp) * self.u0_hat

		return


	def getamp(self):
		tau = (self.t-self.t0)/ self.tE 
		u0 = self.u0.reshape(1, len(self.u0))
		thetaE_hat = self.thetaE_hat.reshape(1,len(self.thetaE_hat))
		tau = tau.reshape(len(tau), 1)

		#shape of u [n, 2]
		u = u0 + tau * thetaE_hat

		#shape of u amp
	#	pdb.set_trace()

		u_amp = np.apply_along_axis(np.linalg.norm, 1, u)
		A = (u_amp**2 + 2)/(u_amp*np.sqrt(u_amp**2 +4))
		return A

def testPSPL():
	#initial input variables 
	imL = 10.0 #solar mass 
	idL = 4000.0 #pc
	idS = 8000.0 #pc
	t0 = 57000.0
	tr = 2000.0
	muS =  np.array([8.0, 0.0])
	muL =  np.array([0.00, 0.00])
	beta = 1.8
	x0S = 0.0
	y0S = 0.0
	y0L = 100.0

	draw_PSPL(imL, idL, idS, t0, tr, muS, muL, beta, x0S, y0S, y0L)

	return

def draw_PSPL(imL, idL, idS, t0, tr, muS, muL, beta, x0S, y0S, y0L):
	ac = PSPL(imL, idL, idS, t0, tr, muS, muL, beta, x0S, y0S, y0L)

	#--------------DISPLAY---------------------------------------
	origin = vector(0,0, -ac.idL)
	#observers
	OBPos = vector(0,0,1)
	OB = sphere(pos = OBPos, radius = 1, color = white)

	#SOURCE positioning
	sPos = vector(ac.x0S, ac.y0S, -ac.idL)+OBPos
	SOURCE = sphere(pos = sPos, radius = 30, color=blue)

	#LENS positioning 
	lPos = vector(ac.x0L, ac.y0L, -ac.idS)+OBPos
	LENS = sphere(pos=lPos, radius = 50, color = yellow)
	LENS.velocity = vector(1, 0, 0)
	ER = ring(pos = lPos, radius = ac.eradiuspc_adjusted, axis = (0,0,1), thickness=15, color= white)

	#LETS DO MATH FOR THE LIGHT CURVES
	thetaS = diff_angle(LENS.pos, SOURCE.pos)
	lthetaplus = (thetaS+np.sqrt((thetaS**2)+(4*ac.thetaE)))/2
	lthetaminus = (thetaS-np.sqrt((thetaS**2)+(4*ac.thetaE)))/2

	thetaS1 = thetaS * radtomas
	lthetaplus1 = lthetaplus* radtomas
	lthetaminus1 = lthetaminus * radtomas

	opacityplus = 0.3
	opacityminus = 0.1
	opacitycen = (lthetaplus*opacityplus+lthetaminus*opacityminus)/(opacityplus+opacityminus)

	ldistplusm = np.tan(lthetaplus)*ac.dLS
	ldistpluspc = ldistplusm * mtopc
	ldistpluspc_adjusted = ldistpluspc/(ac.thetaE1*.8)

	ldistminusm = np.tan(lthetaminus)*ac.dLS
	ldistminuspc = ldistminusm * mtopc
	ldistminuspc_adjusted = ldistminuspc*(ac.thetaE1*10**5)


	cent_adjusted = (ldistpluspc_adjusted*opacityplus+ldistminuspc_adjusted*opacityminus)/(opacityplus+opacityminus)

	#LETS DRAW THE LIGHT CURVES
	plpos = vector(-ldistpluspc_adjusted,0,-ac.idL)
	pluslight = sphere(pos=plpos, radius = 20, color = yellow, opacity = opacityplus)
	plneg = vector(-ldistminuspc_adjusted, 0, -ac.idL)
	minuslight = sphere(pos=plneg, radius = 20, color = red, opacity = opacityminus)
	cenpos = vector(-cent_adjusted, 0, -ac.idL)
	cenlight = sphere(pos=cenpos, radius = 20, color = white, opacity = opacitycen)
	#LABELS
	lmasslabel = label(pos = origin, text = 'ML: '+ str(ac.imL) +' solar masses', xoffset = xoff, yoffset = 140, height = textsize, color = white, line = False)
	ldistancelabel = label(pos = origin, text = 'DL: '+str(ac.idL)+' parsecs', xoffset = xoff,yoffset= 120, height = textsize, color = white, line = False)
	sdistancelabel = label(pos = origin, text = 'DS: '+str(ac.idS)+ ' parsecs', xoffset = xoff, yoffset = 100, height = textsize, color = white, line = False)
	erlabel1 = label(pos = origin, text = 'ER(mas): '+str(ac.thetaE1), xoffset = xoff, yoffset = 80, height = textsize, color =white, line = False)

	llabel = label(pos = LENS.pos, text = 'L', height = textsize-2, color =white, line = False)
	slabel = label(pos = SOURCE.pos, text = 'S', yoffset = 2, height = textsize, color =white, line = False)
	erlabel = label(pos = ER.pos, text = 'ER', yoffset = llabel.yoffset+50, height = textsize-2, color = white, line = False)

	amplabel = label(pos = origin, yoffset = -120, text = '', height = textsize, line=False)
	tlabel = label(pos=origin, yoffset= -100, text = '', height = textsize, line = False)
	sourceposlabel = label(pos = origin, yoffset = -100, xoffset = xoff, text = '', height = textsize, line=False)
	lensposlabel = label(pos = origin, yoffset = -120, xoffset = xoff, text = '', height = textsize, line=False)


	#DISPLAY
	A = ac.getamp()
	lvel = LENS.velocity
	x = 0
	#GRAPHS
	ampgraph = gdisplay(x=0, y = 350, width=500, height=300, title = 'AMP vs T', xtitle = 't', ytitle = 'amp', ymin = 1, ymax = A[ac.tr], xmin = ac.t0-ac.tr, xmax = ac.t0+ac.tr)
	ampcurve = gcurve(gdisplay = ampgraph, color = white)

	lightcoordinatesx = gdisplay(x=800, y = 350, width=500, height=300, title = 'LIGHT CURVE(X) vs T', xtitle = 't', ytitle = 'x value',xmin = ac.t0-ac.tr, xmax = ac.t0+ac.tr)
	pluslightx = gcurve(gdisplay = lightcoordinatesx, color = yellow)
	minuslightx = gcurve(gdisplay = lightcoordinatesx, color = red)
	cenlightx = gcurve(gdisplay = lightcoordinatesx, color = white)

	lightcoordinatesy = gdisplay(x=800, y = 0, width=500, height=300, title = 'LIGHT CURVE(Y) vs T', xtitle = 't', ytitle = 'y value',xmin = ac.t0-ac.tr, xmax = ac.t0+ac.tr)
	pluslighty = gcurve(gdisplay = lightcoordinatesy, color = yellow)
	minuslighty = gcurve(gdisplay = lightcoordinatesy, color = red)
	cenlighty = gcurve(gdisplay = lightcoordinatesy, color = white)
	
	lightxy = gdisplay(x=550, y = 350, width =300, height = 300, title = 'LIGHT CURVE (X VS Y)', xtitle = 'x', ytitle = 'y')
	pluslightxy = gcurve(gdisplay = lightxy, color = yellow)
	minuslightxy = gcurve(gdisplay = lightxy, color = red)
	cenlightxy = gcurve(gdisplay = lightxy, color = white)
	
	for time in ac.t:
		LENS.pos = LENS.pos + lvel
		llabel.pos = llabel.pos + lvel
		ER.pos = ER.pos + lvel
		erlabel.pos = erlabel.pos + lvel
		xdiff = LENS.pos.x-SOURCE.pos.x
		if xdiff !=0:		
			rotateangle = np.arctan(200.0/xdiff)
		else:
			rotateangle = 0.0
		if time-ac.t0 > -ac.tr and time-ac.t0 < 0.0:
			pluslight.pos.x = -np.cos(rotateangle)*ldistpluspc_adjusted
			pluslight.pos.y = -np.sin(rotateangle)*ldistpluspc_adjusted
			minuslight.pos.x = -np.cos(rotateangle)*ldistminuspc_adjusted
			minuslight.pos.y = -np.sin(rotateangle)*ldistminuspc_adjusted
			cenlight.pos.x = -np.cos(rotateangle)*cent_adjusted
			cenlight.pos.y = -np.sin(rotateangle)*cent_adjusted
		elif time-ac.t0 > 0.0 and time-ac.t0 < ac.tr:
			pluslight.pos.x = np.cos(rotateangle)*ldistpluspc_adjusted
			pluslight.pos.y = np.sin(rotateangle)*ldistpluspc_adjusted
			minuslight.pos.x = np.cos(rotateangle)*ldistminuspc_adjusted
			minuslight.pos.y = np.sin(rotateangle)*ldistminuspc_adjusted
			cenlight.pos.x = np.cos(rotateangle)*cent_adjusted
			cenlight.pos.y = np.sin(rotateangle)*cent_adjusted


		pluslight.opacity = opacityplus*A[x]
		minuslight.opacity = opacityminus*A[x]
		cenlight.opacity = opacitycen*A[x]

		amplabel.text = 'AMP: '+str(A[x])+' (scaled e80)'+' REAL AMP: '+str(A[x])
		tlabel.text = 'Time: '+str(time-ac.t0)	
		sourceposlabel.text = 'source x: '+str(SOURCE.pos.x)+' y: '+str(SOURCE.pos.y)
		lensposlabel.text = 'lens x: '+str(LENS.pos.x)+' y: '+str(LENS.pos.y)

		#update graphs
		ampcurve.plot(pos=(time, A[x]))
		pluslightx.plot(pos=(time, pluslight.pos.x))
		pluslighty.plot(pos=(time, pluslight.pos.y))
		pluslightxy.plot(pos=(pluslight.pos.x, pluslight.pos.y))
		minuslightx.plot(pos=(time, minuslight.pos.x))
		minuslighty.plot(pos=(time,minuslight.pos.y))
		minuslightxy.plot(pos=(minuslight.pos.x, minuslight.pos.y))
		cenlightx.plot(pos=(time, cenlight.pos.x))
		cenlighty.plot(pos=(time, cenlight.pos.y))
		cenlightxy.plot(pos=(cenlight.pos.x, cenlight.pos.y))

		x = x+1
		
		print(cent_adjusted)
		rate(500)
		



testPSPL()