
from visual import *
import numpy as np
import pdb
#setup
scene = display(title = 'microlensing', width = 1000, height = 500, background = color.black)

scene.autoscale = False


#random variables
width = -300
height = 200
origin = vector(0,0,0)
textsize = 12

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
opacity = 0.025

#initial input variables 
imL = 5.0 #solar mass 
idL = 4000.0 #pc
idS = 8000.0 #pc
t0 = 57000.0
muS =  np.array([8.0, 0.0])
muL =  np.array([0.00, 0.00])
beta = 1.8
t = np.arange(t0-2000.0, t0+2000.0)
#xS0 = np.arange(-2000,2000)

# conversion
mL = imL * (1.99 * (10 ** 30))
dL = idL * (30856775714409184)
dS = idS * (30856775714409184)
G = 6.67 * 10**-11
c = 3.0 * 10**8
radtomas = 206265000
mtopc = 3.0856776000000000 * 10**16


#relative velocity vector
muRel = muS - muL
muRel1 = np.linalg.norm(muRel)

#calculating einstein radius 
inv_dist_diff = (1.0/dL)-(1.0/dS)
thetaE = (np.sqrt((4.0*G*mL/c**2) * inv_dist_diff)) # radians
thetaE1 = thetaE * radtomas #in mas
thetaE_hat = muRel/ thetaE1

#finding einstein radius distance from center
degreeE = (thetaE * 180)/np.pi # radians to degrees
eradiusm = np.tan(degreeE)*dL #meters
eradiuspc = eradiusm/mtopc
eradiuspc_adjusted = eradiuspc*120000


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



#--------------DISPLAY SHIT---------------------------------------
#observers
OBPos = vector(0,0,1)
OB = sphere(pos = OBPos, radius = 1, color = white)

#lens positioning
lPos = vector(0, 0, -idL)+OBPos
LENS = sphere(pos = lPos, radius = 50, color=blue)
ER = ring(pos = lPos, radius = eradiuspc_adjusted, axis = (0,0,1), thickness=15, color= white)

#source positioning 
sPos = vector(-2000.0, 200.0, -idS)+OBPos
STAR = sphere(pos=sPos, radius = 30, color = yellow)
STAR.velocity = vector(1, 0, 0)

#LETS CREATE THE LIGHT CURVES
plpos = vector(-eradiuspc_adjusted,0,-idL)+OBPos
pluslight = sphere(pos=plpos, radius = 20, color = orange, opacity = opacity)
plneg = vector(eradiuspc_adjusted, 0, -idL)+OBPos
minuslight = sphere(pos=plneg, radius = 20, color = orange, opacity = opacity)

#LABELS
lmasslabel = label(pos = LENS.pos, text = 'lens mass: '+ str(imL) +' solar masses', xoffset = width, yoffset = height, height = textsize, color = white, line = False, align = 'right')
ldistancelabel = label(pos = LENS.pos, text = 'distance to lens: '+str(idL)+' parsecs', xoffset = width+120,yoffset= height, height = textsize, color = white, line = False, align = 'right')
sdistancelabel = label(pos = LENS.pos, text = 'distance to source: '+str(idS)+ ' parsecs', xoffset = width+350, yoffset = height, height = textsize, color = white, line = False, align = 'right')
erlabel1 = label(pos = LENS.pos, text = 'einstein radius (angular): '+str(thetaE1)+ ' MAS', xoffset = width+501, yoffset = height, height = textsize, color =white, line = False, align = 'right')
#erlabel2 = label(pos = LENS.pos, text = '(distance): '+str(eradiuspc)+ ' parsecs', xoffset = width+700, yoffset = height, height = textsize, color = white, line = False, align = 'right')
llabel = label(pos = LENS.pos, text = 'LENS', yoffset=5, height = textsize, color =white, line = False)
slabel = label(pos = STAR.pos, text = 'SOURCE', yoffset=5, height = textsize, color =white, line = False)
erlabel = label(pos = ER.pos, text = 'EINSTEIN RADIUS', yoffset = slabel.yoffset+60, height = textsize-2, color = white, line = False)

amplabel = label(pos = LENS.pos, yoffset = -200, text = 'AMP', height = 10, line=False)
tlabel = label(pos=LENS.pos, yoffset= -100, text = 'time', height = 10, line = False)


#moving the lems

def movingsource():
	deltaT = 0.05
	svel = STAR.velocity
	t = 0
	while t < 20:
		STAR.pos = STAR.pos + svel
		slabel.pos = slabel.pos + svel

		u = diff_angle(STAR.pos, LENS.pos)
		us = u/thetaE
		amp = (us**2+2)/(us*np.sqrt(us**2+4))
		loga = np.log(amp)
			#mas
		u1 = (u *radtomas)
		u1s = u1/thetaE1
		amp1 = ((u1s**2)+2)/(u1s*(np.sqrt((u1s**2)+4)))
		loga1 = np.log(amp1)

		amplabel.text = 'AMP: '+str(loga)
		tlabel.text = 'Time: '+str(t)
		pluslight.opacity = opacity*loga
		minuslight.opacity = opacity*loga

		print(amp1)

		t = t+deltaT
		rate(100)

	
def movingsource1(t = t, t0=t0, A = getamp()):
	svel = STAR.velocity
	x = 0
	A = A**10
	for time in t:
		STAR.pos = STAR.pos + svel
		slabel.pos = slabel.pos + svel

		amplabel.text = 'AMP: '+str(A[x])
		tlabel.text = 'Time: '+str(time-t0)
		pluslight.opacity = opacity*A[x]
		minuslight.opacity = opacity*A[x]

	##	print('%.15f' %amp)
	##	print('%.15f' %amp1)
		x = x+1
		rate(100)


movingsource1()
'''
			#radians
			#BE ABLE TO USE THIS ANGLE TO CREATE A VECTOR USE TRIG????
		thetaS = diff_angle(STAR.pos, LENS.pos) 
		thetaPOS = (thetaS + (np.sqrt((thetaS**2)+(4*thetaE))))/2
		thetaNEG = (thetaS - (np.sqrt((thetaS**2)+(4*thetaE))))/2

			#mas
		thetaS1 = thetaS * radtomas
		thetaPOS1 = (thetaS1 + np.sqrt((thetaS1**2)+(4*thetaE1)))/2
		thetaNEG1= (thetaS1 - np.sqrt((thetaS1**2)+(4*thetaE1)))/2

	'''	