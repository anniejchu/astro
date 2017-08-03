from visual import *
import numpy as np 

scene.width = 600
scene.height = 600
scene.autoscale = False
scene.range = 100
scene.title = 'BLIP'

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
	def __init__(self, imL, idL, idS, t0, muS, muL, beta, x0S, y0S, x0L, y0L):
		self.imL = imL
		self.idL = idL
		self.idS = idS
		self.t0 = t0
		self.muS = muS
		self.muL = muL
		self.beta = beta
		self.x0S = x0S
		self.y0S = x0S
		self.x0L = x0L
		self.y0L = y0L

		#basic
		self.tr = self.x0L-self.x0S
		self.t = np.arange(self.t0-self.tr, self.t0+self.tr)
		self.mL = self.imL * smtokg
		self.dL = self.idL * pctom
		self.dS = self.idS * pctom
		self.dLS = self.dS-self.dL
		self.idLS = self.idS-self.idL
		self.off = (tan(self.beta/radtomas)*self.idLS)*10**5 
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
		self.eradiuspc_adjusted = eradiuspc*10**6.5


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

		#shape of u amp
	#	pdb.set_trace()

		u_amp = np.apply_along_axis(np.linalg.norm, 1, u)
		A = (u_amp**2 + 2)/(u_amp*np.sqrt(u_amp**2 +4))
		return A

	def get_centroid_shift(self):
		dt_in_years = (self.t - self.t0)/365
		tau = (self.t - self.t0) / self.tE

		# Shape of arrays:
		# thetaS: [N_times, 2]
		# u: [N_times, 2]
		# u_amp: [N_times]
		thetaS = self.thetaS0 + np.outer(dt_in_years, self.muRel) # mas
		u = thetaS / self.thetaE1
		u_amp = np.apply_along_axis(np.linalg.norm, 1, u)

		shift_norm_factor = u_amp**2 + 2.0

		shift = thetaS
		shift[:, 0] /= shift_norm_factor
		shift[:, 1] /= shift_norm_factor
                    
		return shift


def randostars():
	rx = []
	ry = []
	opacity = []
	radius=  []
	for x in range(0, 7):
		rx.append(random.randint(-90,90))
		ry.append(random.randint(-90,90))
		opacity.append(random.uniform(0.2, 0.8))
		radius.append(random.uniform(1,2))
		sphere(pos=(rx[x],ry[x],0), radius = radius[x], color = color.white, opacity = opacity[x])


def sourcemoving(imL, idL, idS, t0, muS, muL, beta, x0S, y0S, x0L, y0L):
	it = PSPL(imL, idL, idS, t0, muS, muL, beta, x0S, y0S, x0L, y0L)
	opacity = 0.1

	#####################################################################
	source = sphere(pos=(x0S,y0S,0), radius = 1.5, color = color.white, opacity = 0.5)
	lens = ring(pos=(x0L, y0L, 0), radius = 2.5, thickness = 0.4, color = color.white, axis = (0,0,1))
	centroid = sphere(pos = lens.pos, radius = 1, color = color.white, opacity = opacity)
	sourcev = vector(x0S,y0S,0)

	movex = ((x0L-x0S)-it.off)/(x0L-x0S)
	movey = (y0L-y0S)/(x0L-x0S)
	
	A = it.getamp()**10
	cs = it.get_centroid_shift()*5
	csx = cs[:,0]
	csy = cs[:,1]
	m = 0
	o = len(A)/2-5
	for time in it.t:
		source.x = movex + source.x
		source.y = movey + source.y
		if time-it.t0 < 2*it.off and time-it.t0 > -2*it.off:
			centroid.x = lens.x-csx[m]
			centroid.y = lens.y-csy[m]
			centroid.opacity = opacity * A[o]
			source.opacity = 0.0
			m = m+10
			o = o+5
		elif time-it.t0 > -it.off:
			centroid.opacity = 0.0
			source.opacity = 0.5
	

		rate(10)

def testing():
	imL = 10.0 #solar mass 
	idL = 4000.0 #kpc
	idS = 8000.0 #kpc
	t0 = 57000.0
	muS =  np.array([0.0, 0.0])
	muL =  np.array([8.00, 0.00])
	beta = 1.8
	x0S = -90.0
	y0S = -90.0
	x0L = np.random.randint(-30,30)
	y0L = np.random.randint(-30,30)

	sourcemoving(imL, idL, idS, t0, muS, muL, beta, x0S, y0S, x0L, y0L)

randostars()
testing()