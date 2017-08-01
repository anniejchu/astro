
from visual import *
from visual.graph import *
import numpy as np
import pdb
from astropy import constants as const
from astropy import units
from astropy.time import Time

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
days_per_year = 365.25 

#---old
'''
# conversion
mL = imL * (1.99 * (10 ** 30))
dL = idL * (30856775714409184)
dS = idS * (30856775714409184)
dLS = dS-dL
G = 6.67 * 10**-11
c = 3.0 * 10**8
radtomas = 206265000
mtopc = 1/(3.0856776000000000 * 10**16)

#calculating einstein radius 
inv_dist_diff = (1.0/dL)-(1.0/dS)
thetaE = (np.sqrt((4.0*G*mL/c**2) * inv_dist_diff)) # radians
thetaE1 = thetaE * radtomas #in mas
thetaE_hat = muRel/ thetaE1

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

	u_amp = np.apply_along_axis(np.linalg.norm, 1, u)
	A = (u_amp**2 + 2)/(u_amp*np.sqrt(u_amp**2 +4))
	return A
'''

class PSPL(object): 

	def __init__(self, mL, dL, dS, t0, tr, xS0, beta, muL, muS):
		"""
		INPUTS: 
		mL = mass of lens (solar mass)
		dL = distance from observer to lens (pc)
		ds = distance from observer to source (pc)
		t0 = time of photometric peak (seen from earth)
		xS0 = initial position of the source
		tr = time range from t0 
		beta = angular distance between lens and source on the plane of the skay (mas)
		muL = vector [RA, Dec] lens proper motion (mas/yr)
		muS = vector [RA, Dec] source proper motion (mas/yr)

		"""
		self.mL = mL 
		self.dL = dL
		self.dS = dS
		self.t0 = t0
		self.tr = tr
		self.x0S = x0S
		self.beta = beta
		self.muL = muL
		self.muS = muS

		#calculating the relative parallax
		inv_dist_diff = (1.0 / (dL * units.pc)) - (1.0 / (dS * units.pc))
		piRel = units.rad * units.au * inv_dist_diff
		self.piRel = piRel.to('mas').value

		# Calculate the individual parallax ---NOT USED----
		piS = (1.0 / self.dS) * (units.rad * units.au / units.pc)
		piL = (1.0 / self.dL) * (units.rad * units.au / units.pc)
		self.piS = piS.to('mas').value
		self.piL = piL.to('mas').value     
           
		#relative velocity vector
		self.muRel = self.muS - self.muL
		self.muRel_amp = np.linalg.norm(self.muRel)

		#calculating the einstein radius 
		self.ithetaE = units.rad * np.sqrt((4.0 * const.G * mL * units.M_sun / const.c**2) * inv_dist_diff) #IN RADIANS used in my drawings!!!!!!!!!!!!!!
		self.thetaE_amp = thetaE.to('mas').value #mas
		self.thetaE_hat = self.muRel / self.muRel_amp
		self.thetaE = self.thetaE_amp * self.thetaE_hat #NOT USED

		# Calculate the microlensing parallax (NOT USED)
		self.piE = (self.piRel / self.thetaE_amp) * self.thetaE_hat

        #calculating the closest approach vector
		self.u0_hat = np.zeros(2, dtype=float)
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
                

		self.u0_amp = self.beta / self.thetaE_amp  # in Einstein units
		self.u0 = np.abs(self.u0_amp) * self.u0_hat

		''' ---------------------NOT USED-----------------------
         # Angular separation vector between source and lens (vector from lens to source)
        self.thetaS0 = self.u0 * self.thetaE_amp    # mas

        # Calculate the position of the lens on the sky at time, t0
        self.xL0 = self.xS0 - (self.thetaS0 * 1e-3)

        # Calculate the microlensing parallax
        self.piE_amp = self.piRel / self.thetaE_amp
        self.piE = self.piE_amp * self.thetaE_hat
		--------------------------------------------------------
		'''

        #calculating the einstein crossing time
		self.tE = (self.thetaE_amp / self.muRel_amp) * days_per_year
		return

	def getamp(self):
		self.t = np.arange(self.t0-self.tr, self.t0+self.tr)
		tau = (self.t-self.t0) / self.tE

    	#convert to matrices
    	#	matrix shapes:
    	#		u0, thetaE_hat = [1,2]
    	#		tau = [n times, 1]
		u0 = self.u0.reshape(1, len(self.u0))
		thetaE_hat = self.thetaE_hat.reshape(1, len(self.thetaE_hat))
		tau = tau.reshape(len(tau), 1)

        # Shape of u: [N_times, 2]        
		u = u0 + tau * thetaE_hat

        # Shape of u_amp: [N_times]
		u_amp = np.linalg.norm(u, axis=1)
		A = (u_amp**2 + 2) / (u_amp * np.sqrt(u_amp**2 + 4))
		return A
	'''
    #-----------------NOT USED ---------------------------
  
	def get_centroid_shift(self, t):
        """Get the centroid shift (in mas) for a list of
        observation times (in MJD).
        """
        
        dt_in_years = (t - self.t0) / days_per_year
        tau = (t - self.t0) / self.tE

        # Shape of arrays:
        # thetaS: [N_times, 2]
        # u: [N_times, 2]
        # u_amp: [N_times]
        thetaS = self.thetaS0 + np.outer(dt_in_years, self.muRel) # mas
        u = thetaS / self.thetaE_amp
        u_amp = np.linalg.norm(u, axis=1)

        shift_norm_factor = u_amp**2 + 2.0

        shift = thetaS
        shift[:, 0] /= shift_norm_factor
        shift[:, 1] /= shift_norm_factor
                    
        return shift
	'''

	def einsteinringdis(self):
		#finding einstein radius distance from center #need to add to here code
		eradiusm = np.tan(ithetaE)*(dL* units.pc) #meters
		eradiuspc = eradiusm*(units.pc.to(units.m))
		eradiuspc_adjusted = eradiuspc*10**7

		return eradiuspc_adjusted


#--------------DISPLAY---------------------------------------
'''
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
ldistminuspc_adjusted = ldistminuspc*(thetaE1*10**5)


cent_adjusted = (ldistpluspc_adjusted+ldistminuspc_adjusted)/2

#LETS CREATE THE LIGHT CURVES
plpos = vector(-ldistpluspc_adjusted,0,-idL)+OBPos
pluslight = sphere(pos=plpos, radius = 20, color = cyan, opacity = opacityplus)
plneg = vector(-ldistminuspc_adjusted, 0, -idL)+OBPos
minuslight = sphere(pos=plneg, radius = 20, color = red, opacity = opacityminus)
cenpos = vector(-cent_adjusted, 0, -idL)
cenlight = sphere(pos=cenpos, radius = 20, color = white, opacity = opacitycen)
#LABELS
lmasslabel = label(pos = origin, text = 'ML: '+ str(imL) +' solar masses', xoffset = xoff, yoffset = 140, height = textsize, color = white, line = False)
ldistancelabel = label(pos = origin, text = 'DL: '+str(idL)+' parsecs', xoffset = xoff,yoffset= 120, height = textsize, color = white, line = False)
sdistancelabel = label(pos = origin, text = 'DS: '+str(idS)+ ' parsecs', xoffset = xoff, yoffset = 100, height = textsize, color = white, line = False)
erlabel1 = label(pos = origin, text = 'ER(mas): '+str(thetaE1), xoffset = xoff, yoffset = 80, height = textsize, color =white, line = False)

llabel = label(pos = LENS.pos, text = 'L', height = textsize-2, color =white, line = False)
slabel = label(pos = SOURCE.pos, text = 'S', yoffset = 2, height = textsize, color =white, line = False)
erlabel = label(pos = ER.pos, text = 'ER', yoffset = llabel.yoffset+50, height = textsize-2, color = white, line = False)

amplabel = label(pos = origin, yoffset = -120, text = '', height = textsize, line=False)
tlabel = label(pos=origin, yoffset= -100, text = '', height = textsize, line = False)
sourceposlabel = label(pos = origin, yoffset = -100, xoffset = xoff, text = '', height = textsize, line=False)
lensposlabel = label(pos = origin, yoffset = -120, xoffset = xoff, text = '', height = textsize, line=False)



def moving(t = t, t0=t0, A = getamp()):
	lvel = LENS.velocity
	x = 0
	A1 = A**5
	#GRAPHS
	ampgraph = gdisplay(x=0, y = 350, width=500, height=300, title = 'AMP vs T', xtitle = 't', ytitle = 'amp', ymin = 1, ymax = A1[tr], xmin = t0-tr, xmax = t0+tr)
	ampcurve = gcurve(gdisplay = ampgraph, color = white)
	ampadjcurve = gcurve(gdisplay = ampgraph, color = white)

	lightcoordinatesx = gdisplay(x=800, y = 350, width=500, height=300, title = 'LIGHT CURVE(X) vs T', xtitle = 't', ytitle = 'x value',xmin = t0-tr, xmax = t0+tr)
	pluslightx = gcurve(gdisplay = lightcoordinatesx, color = cyan)
	minuslightx = gcurve(gdisplay = lightcoordinatesx, color = red)
	cenlightx = gcurve(gdisplay = lightcoordinatesx, color = white)

	lightcoordinatesy = gdisplay(x=800, y = 0, width=500, height=300, title = 'LIGHT CURVE(Y) vs T', xtitle = 't', ytitle = 'y value',xmin = t0-tr, xmax = t0+tr)
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
		


moving()
'''
print(ithetaE)