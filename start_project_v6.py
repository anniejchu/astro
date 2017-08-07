from visual import *
from visual.graph import *
import numpy as np
import pdb
#setup
scene = display(title = 'microlensing', width = 750, height = 350, background = color.black)
scene.autoscale = False
scene.range = 600

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
	def __init__(self, imL, idL, idS, t0, tr, muS, muL, beta, x0S, y0S):
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

		#basic
		self.t = np.arange(self.t0-self.tr, self.t0+self.tr)
		self.mL = self.imL * smtokg
		self.dL = self.idL * pctom
		self.dS = self.idS * pctom
		self.idLS = self.idS-self.idL
		self.dLS = self.dS-self.dL
		self.x0L = -self.tr
		self.y0L = (tan(self.beta/radtomas)*self.idLS)*100
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
		self.eradiuspc = eradiusm*mtopc #pc
		#self.eradiuspc_adjusted = self.eradiuspc*10**6.8
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

		#angular separation between source and lens vector
		self.thetaS0 = self.u0 * self.thetaE1
		
		return

	def getamp(self):
		tau = (self.t-self.t0)/ self.tE 
		u0 = self.u0.reshape(1, len(self.u0))
		thetaE_hat = self.thetaE_hat.reshape(1,len(self.thetaE_hat))
		tau = tau.reshape(len(tau), 1)
		#shape of u [n, 2]
		u = u0 + tau * thetaE_hat
		u_amp = np.apply_along_axis(np.linalg.norm, 1, u)
		A = (u_amp**2 + 2)/(u_amp*np.sqrt(u_amp**2 +4))
		
		return A

	def get_centroid_shift(self):
		dt_in_years = (self.t - self.t0)/365
		tau = (self.t - self.t0) / self.tE

		# Shape of arrays: 1) thetaS: [N_times, 2] 2) u: [N_times, 2] 3) u_amp: [N_times]
		thetaS = self.thetaS0 + np.outer(dt_in_years, self.muRel) # mas
		u = thetaS / self.thetaE1
		u_amp = np.apply_along_axis(np.linalg.norm, 1, u)

		shift_norm_factor = u_amp**2 + 2.0

		shift = thetaS
		shift[:, 0] /= shift_norm_factor
		shift[:, 1] /= shift_norm_factor                    
		
		return shift

def testPSPL():
	#initial input variables 
	imL = float(input("Enter BH mass (in solar masses): ")) #solar mass 5.0
	idL = float(input("Enter distance to lens (pc): ")) #4000.0 #pc
	idS = float(input("Enter distance to source (pc): "))#8000.0 #pc
	t0 = float(input("enter t0: ")) #57000.0
	tr = float(input("enter time range: ")) #2000.0
	muS =  np.array([0.00, 0.00])
	muL =  np.array([8.00, 0.00])
	beta = float(input("enter beta: ")) #2.0
	x0S = 0.0
	y0S = 0.0

	draw_PSPL(imL, idL, idS, t0, tr, muS, muL, beta, x0S, y0S)

	return



def draw_PSPL(imL, idL, idS, t0, tr, muS, muL, beta, x0S, y0S):
	ac = PSPL(imL, idL, idS, t0, tr, muS, muL, beta, x0S, y0S)
	rA = ac.getamp()
	rcs = ac.get_centroid_shift()
	A = ac.getamp()**10
	cs = ac.get_centroid_shift()*100
	csx = cs[:, 0]
	csy = cs[:, 1]
	plc = ac.get_centroid_shift()*650
	plcx = plc[:,0]
	plcy = plc[:,1]
	mlc = -ac.get_centroid_shift()*200
	mlcx = mlc[:,0]
	mlcy = mlc[:,1]
	eradiuspc_adjusted = ac.eradiuspc*(min(plcy)*10000)
	print(ac.eradiuspc)
	#--------------DISPLAY---------------------------------------
	origin = vector(0,0, -ac.idL)
	#observers

	#SOURCE positioning
	sPos = vector(ac.x0S, ac.y0S, -ac.idS)
	SOURCE = sphere(pos = sPos, radius = 30, color=blue)

	#LENS positioning 
	lPos = vector(ac.x0L, ac.y0L, -ac.idL)
	LENS = sphere(pos=lPos, radius = 20, color = yellow)
	LENS.velocity = vector(1, 0, 0)

	erpos = vector(lPos.x, lPos.y, -ac.idL)
	ER = ring(pos = erpos, radius = eradiuspc_adjusted, axis = (0,0,1), thickness=15, color= white, opacity = 0.5)
	
	#LETS DO MATH FOR THE LIGHT CURVES
	thetaS = diff_angle(LENS.pos, SOURCE.pos)
	lthetaplus = (thetaS+np.sqrt((thetaS**2)+(4*ac.thetaE)))/2
	lthetaminus = (thetaS-np.sqrt((thetaS**2)+(4*ac.thetaE)))/2

	thetaS1 = thetaS * radtomas
	lthetaplus1 = lthetaplus* radtomas
	lthetaminus1 = lthetaminus * radtomas

	opacityplus = 0.03
	opacityminus = 0.01
	opacitycen = (lthetaplus*opacityplus+lthetaminus*opacityminus)/(lthetaplus+lthetaminus)

	#LETS DRAW THE LIGHT CURVES
	plpos = vector(SOURCE.pos)
	pluslight = sphere(pos=plpos, radius = 30, color = yellow, opacity = opacityplus)
	plneg = vector(SOURCE.pos)
	minuslight = sphere(pos=plneg, radius = 30, color = red, opacity = opacityminus)
	cenpos = vector(SOURCE.pos)
	cenlight = sphere(pos=cenpos, radius = 30, color = white, opacity = opacitycen)
	#LABELS
	plabel = label(pos = origin, text = 'PARAMETERS', xoffset = xoff-5, yoffset = 140, height = textsize+3, color = white, line = False)
	lmasslabel = label(pos = origin, text = 'ML: '+ str(ac.imL) +' solar masses', xoffset = xoff, yoffset = 110, height = textsize, color = white, line = False, box = False)
	ldistancelabel = label(pos = origin, text = 'DL: '+str(ac.idL)+' parsecs', xoffset = xoff,yoffset= 90, height = textsize, color = white, line = False, box = False)
	sdistancelabel = label(pos = origin, text = 'DS: '+str(ac.idS)+ ' parsecs', xoffset = xoff, yoffset = 70, height = textsize, color = white, line = False, box = False)
	betalabel = label(pos = origin, text = 'beta: '+str(ac.beta)+ ' mas', xoffset = xoff, yoffset = 50, height = textsize, color = white, line = False, box = False)
	tr0label = label(pos = origin, text = 't0: '+str(ac.t0)+' tr: '+str(ac.tr), xoffset = xoff, yoffset = 30, height = textsize, color = white, line = False, box = False)


	line = paths.line(start=(-550, 30), end=(-300, 30))
	curve(pos=line.pos)

	mlalabel = label(pos = origin, text = 'MICROLENSING ANIMATOR', xoffset=160, yoffset = 140, height = textsize+3, line=False)

	erlabel1 = label(pos = origin, text = 'ER(mas): '+str(round(ac.thetaE1, 2))+' mas', xoffset = xoff, yoffset = 0, height = textsize, color =white, line = False, box = False)
	telabel= label(pos = origin, text = 'TE: '+str(int(ac.tE))+' days', xoffset = xoff, yoffset= -20, height = textsize, color = color.white, line = False, box = False)
	'''
	llabel = label(pos = LENS.pos, text = 'L', height = textsize-2, color =white, line = False)
	slabel = label(pos = SOURCE.pos, text = 'S', yoffset = 2, height = textsize, color =white, line = False)
	erlabel = label(pos = ER.pos, text = 'ER', yoffset = llabel.yoffset+50, height = textsize-2, color = white, line = False)
	'''
	amplabel = label(pos = origin, yoffset = -130, text = '', height = textsize+2, line=False, box = False)
	tlabel = label(pos=origin, yoffset= -110, text = '', height = textsize+2, line = False, box = False)
	sourceposlabel = label(pos = origin, yoffset = -110, xoffset = -xoff, text = '', height = textsize, line=False, box = False)
	lensposlabel = label(pos = origin, yoffset = -130, xoffset = -xoff, text = '', height = textsize, line=False, box = False)


	#DISPLAY

	lvel = LENS.velocity
	x = 0
	#GRAPHS
	ampgraph = gdisplay(x=0, y = 340, width=500, height=500, title = 'AMP vs T', xtitle = 't', ytitle = 'amp', ymin = 1, ymax = rA[ac.tr], xmin = ac.t0-ac.tr, xmax = ac.t0+ac.tr)
	ampcurve = gcurve(gdisplay = ampgraph, color = white)

	lightcoordinatesx = gdisplay(x=880, y = 340, width=500, height=500, title = 'LIGHT CURVE(X) vs T', xtitle = 't', ytitle = 'x value', ymin = min(plcx), ymax = max(plcx), xmin = ac.t0-ac.tr, xmax = ac.t0+ac.tr)
	pluslightx = gcurve(gdisplay = lightcoordinatesx, color = yellow)
	minuslightx = gcurve(gdisplay = lightcoordinatesx, color = red)
	cenlightx = gcurve(gdisplay = lightcoordinatesx, color = white)

	lightcoordinatesy = gdisplay(x=725, y = 0, width=650, height=350, title = 'LIGHT CURVE(Y) vs T', xtitle = 't', ytitle = 'y value', ymin = plcy[ac.tr], ymax = -plcy[ac.tr], xmin = ac.t0-ac.tr, xmax = ac.t0+ac.tr)
	pluslighty = gcurve(gdisplay = lightcoordinatesy, color = yellow)
	minuslighty = gcurve(gdisplay = lightcoordinatesy, color = red)
	cenlighty = gcurve(gdisplay = lightcoordinatesy, color = white)
	
	lightxy = gdisplay(x=480, y = 340, width =420, height = 500, title = 'LIGHT CURVE (X VS Y)', xtitle = 'x', ytitle = 'y', ymin = min(plcy), ymax = -min(plcy), xmin = min(plcx), xmax = max(plcx))
	pluslightxy = gcurve(gdisplay = lightxy, color = yellow)
	minuslightxy = gcurve(gdisplay = lightxy, color = red)	
	cenlightxy = gcurve(gdisplay = lightxy, color = white)
	
	for time in ac.t:
		LENS.pos.x = LENS.pos.x + lvel.x
		#llabel.pos = llabel.pos + lvel
		ER.pos.x = ER.pos.x + lvel.x
		#erlabel.pos = erlabel.pos + lvel
		
		cenlight.pos.x = csx[x]
		cenlight.pos.y = csy[x]
		pluslight.pos.x = plcx[x]
		pluslight.pos.y = plcy[x]
		minuslight.pos.x = mlcx[x]
		minuslight.pos.y = mlcy[x]


		pluslight.opacity = opacityplus*A[x]
		minuslight.opacity = opacityminus*A[x]
		cenlight.opacity = opacitycen*A[x]
		amplabel.text = 'AMP : '+str(round(rA[x], 6))
		tlabel.text = 'Time (days): '+str(time)	
		sourceposlabel.text = 'source x: '+str(SOURCE.pos.x)+' y: '+str(SOURCE.pos.y)
		lensposlabel.text = 'lens x: '+str(round(LENS.pos.x, 4))+' y: '+str(round(LENS.pos.y, 6))

		#update graphs
		ampcurve.plot(pos=(time, rA[x]))
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
		
		rate(200)
		


testPSPL()