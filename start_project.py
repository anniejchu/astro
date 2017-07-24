
from visual import *
import numpy as np
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
opacity = 0.5

#initial input variables 
imL = 5 #solar mass 
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


#observers
OBPos = vector(0,0,1)
OB = sphere(pos = OBPos, radius = 1, color = white)

#lens positioning
lPos = vector(0, 0, -idL)+OBPos
LENS = sphere(pos = lPos, radius = 50, color=blue)
ER = ring(pos = lPos, radius = eradiuspc_adjusted, axis = (0,0,1), thickness=15, color= white)

#source positioning 
sPos = vector(-3000.0, 200.0, -idS)+OBPos
STAR = sphere(pos=sPos, radius = 30, color = yellow)
STAR.velocity = vector(250, 0, 0)

#LETS CREATE THE LIGHT CURVES
plpos = vector(-eradiuspc_adjusted,0,-idL)+OBPos
pluslight = sphere(pos=plpos, radius = 20, color = orange, opacity = opacity)
plneg = vector(eradiuspc_adjusted, 0, -idL)+OBPos
minuslight = sphere(pos=plneg, radius = 20, color = orange, opacity = opacity)



#--------------DISPLAY SHIT---------------------------------------

lmasslabel = label(pos = LENS.pos, text = 'lens mass: '+ str(imL) +' solar masses', xoffset = width, yoffset = height, height = textsize, color = white, line = False, align = 'right')
ldistancelabel = label(pos = LENS.pos, text = 'distance to lens: '+str(idL)+' parsecs', xoffset = width+120,yoffset= height, height = textsize, color = white, line = False, align = 'right')
sdistancelabel = label(pos = LENS.pos, text = 'distance to source: '+str(idS)+ ' parsecs', xoffset = width+350, yoffset = height, height = textsize, color = white, line = False, align = 'right')
erlabel1 = label(pos = LENS.pos, text = 'einstein radius (angular): '+str(thetaE1)+ ' MAS', xoffset = width+501, yoffset = height, height = textsize, color =white, line = False, align = 'right')
#erlabel2 = label(pos = LENS.pos, text = '(distance): '+str(eradiuspc)+ ' parsecs', xoffset = width+700, yoffset = height, height = textsize, color = white, line = False, align = 'right')
llabel = label(pos = LENS.pos, text = 'LENS', yoffset=5, height = textsize, color =white, line = False)
slabel = label(pos = STAR.pos, text = 'SOURCE', yoffset=5, height = textsize, color =white, line = False)
erlabel = label(pos = ER.pos, text = 'EINSTEIN RADIUS', yoffset = slabel.yoffset+60, height = textsize-2, color = white, line = False)

'''

#light curve info ASK NIJAID WHAT THIS ANGLE IS
	#radians
	#BE ABLE TO USE THIS ANGLE TO CREATE A VECTOR USE TRIG????
thetaS = diff_angle(STAR.pos, OB.pos) 
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
'''

#moving the lems
def movingsource(t, deltaT, timelimit):
	#FIND EINSTEIN CROSSING TIME, MAKE IT RELEVANT TO THIS FUNCTION
	#REMEMBER THAT THE EINSTEIN CROSSING TIME IS JUST HALF OF THE GRAPH
	while t < timelimit:
		STAR.pos = STAR.pos + STAR.velocity*deltaT
		slabel.pos = slabel.pos + STAR.velocity*deltaT
		#light curve info ASK NIJAID WHAT THIS ANGLE IS
			#radians
			#BE ABLE TO USE THIS ANGLE TO CREATE A VECTOR USE TRIG????
		thetaS = diff_angle(STAR.pos, OB.pos) 
		thetaPOS = (thetaS + np.sqrt((thetaS**2)+(4*thetaE)))/2
		thetaNEG = (thetaS - np.sqrt((thetaS**2)+(4*thetaE)))/2

			#mas
		thetaS1 = thetaS * radtomas
		thetaPOS1 = (thetaS1 + np.sqrt((thetaS1**2)+(4*thetaE1)))/2
		thetaNEG1= (thetaS1 - np.sqrt((thetaS1**2)+(4*thetaE1)))/2


		#light magnification 
  			#radians
		u = diff_angle(STAR.pos, LENS.pos)
		amp = ((u**2)+2)/(u*(np.sqrt((u**2)+4)))
			#mas
		u1 = u *radtomas
		amp1 = ((u1**2)+2)/(u1*(np.sqrt((u1**2)+4)))

		t = t+deltaT
		rate(100)








movingsource(0, 0.05, 20)