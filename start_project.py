
from visual import *
import numpy as np
#setup
scene.width = 1000
scene.height = 500
scene.autoscale = False


#random variables
width = -40
height = 180
origin = (0,0,0)
textsize = 12


#variables 
imL = 5 #solar mass MUST CONVERT
idL = 4000 #pc
idS = 8000 #pc

# conversion
mL = imL * (1.99 * (10 ** 30))
dL = idL * (30856775714409184)
dS = idS * (30856775714409184)
G = 6.67 * 10**-11
c = 3.0 * 10**8
radtomas = 206265000
mtopc = 3.0856776000000000 * 10**16


#calculating einstein radius 
inv_dist_diff = (1.0/dL)-(1.0/dS)
thetaE = (np.sqrt((4.0*G*mL/c**2) * inv_dist_diff)) # radians
thetaE1 = thetaE * radtomas #in mas

#finding einstein radius distance from center
degreeE = (thetaE * 180)/np.pi # radians to degrees
eradiusm = np.tan(degreeE)*dL #meters
eradiuspc = eradiusm/mtopc
eradiuspc_adjusted = eradiuspc*120000

#lens positioning
lPos = vector(0, 0, -idL)
LENS = sphere(pos = lPos, radius = 50, color=color.blue)
ER = ring(pos = lPos, radius = eradiuspc_adjusted, axis = (0,0,1), thickness=15, color= color.white)

#source positioning 
sPos = vector(-3000.0, 200.0, -idS)
STAR = sphere(pos=sPos, radius = 50, color = color.yellow)
STAR.velocity = vector(250, 0, 0)


#deltaT = speed , xlimit = 4000 
'''
def movingsource(deltaT, xlimit):
	while STAR.pos.x < xlimit: 
		STAR.pos = STAR.pos + STAR.velocity*deltaT
		rate(100)
'''
def movingsource(t, deltaT, timelimit):
	#FIND EINSTEIN CROSSING TIME, MAKE IT RELEVANT TO THIS FUNCTION
	#REMEMBER THAT THE EINSTEIN CROSSING TIME IS JUST HALF OF THE GRAPH
	while t < timelimit:
		STAR.pos = STAR.pos + STAR.velocity*deltaT
		slabel.pos = slabel.pos + STAR.velocity*deltaT
		t = t+deltaT
		rate(100)

#--------------DISPLAY SHIT---------------------------------------

lmasslabel = label(pixel_pos = vector(0,20,0), text = 'lens mass: '+ str(imL) +' solar masses', height = textsize, color = color.white, line = False, align = 'left')
ldistancelabel = label(pos = STAR.pos, text = 'distance to lens: '+str(idL)+' parsecs', xoffset = width,yoffset= height-25, height = textsize, color = color.white, line = False, align = 'left')
sdistancelabel = label(pos = STAR.pos, text = 'distance to source: '+str(idS)+ ' parsecs', xoffset = width, yoffset = height-50, height = textsize, color = color.white, line = False, align = 'left')
erlabel1 = label(pos = STAR.pos, text = 'einstein radius (angular): '+str(thetaE1)+ ' MAS', xoffset = width, yoffset = height-325, height = textsize, color = color.white, line = False, align = 'left')
erlabel2 = label(pos = STAR.pos, text = '(distance): '+str(eradiuspc)+ ' parsecs', xoffset = width, yoffset = height-350, height = textsize, color = color.white, line = False, align = 'left')
llabel = label(pos = LENS.pos, text = 'LENS', yoffset=5, height = textsize, color = color.white, line = False)
slabel = label(pos = STAR.pos, text = 'SOURCE', yoffset=5, height = textsize, color = color.white, line = False)
erlabel = label(pos = ER.pos, text = 'EINSTEIN RADIUS', yoffset = slabel.yoffset+60, height = textsize-2, color = color.white, line = False)




#light curve info ASK NIJAID WHAT THIS ANGLE IS
	#radians
	#BE ABLE TO USE THIS ANGLE TO CREATE A VECTOR USE TRIG????
thetaS = diff_angle(STAR.pos, origin) 
thetaPOS = (thetaS + np.sqrt((thetaS**2)+(4*thetaE)))/2
thetaNEG = (thetaS - np.sqrt((thetaS**2)+(4*thetaE)))/2

	#mas
thetaS1 = thetaS * radtomas
thetaPOS1 = (thetaS1 + np.sqrt((thetaS1**2)+(4*thetaE1)))/2
thetaNEG1= (thetaS1 - np.sqrt((thetaS1**2)+(4*thetaE1)))/2


#light magnification 
#USE THIS TO EITHER CHANGE THE COLOR OR INTENSITY(OPAQUENESS)
	#radians
u = diff_angle(STAR.pos, LENS.pos)
amp = ((u**2)+2)/(u*(np.sqrt((u**2)+4)))
	#mas
u1 = u *radtomas
amp1 = ((u1**2)+2)/(u1*(np.sqrt((u1**2)+4)))



movingsource(-20, 0.05, 50)