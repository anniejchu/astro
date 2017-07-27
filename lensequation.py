from visual import *
import numpy as np

ox = -80.0
oy = 0.0
oz = -40.0
opos = vector(ox, oy, oz)

lx = 0.0
ly = 0.0
lz = -40.0
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

#IMAGE STUFF
ok = sz-oz
theta = atan2(ok,dl)
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
observerlabel = label(pos = opos, text='observer', yoffset = -2, height = 10, color = color.white, line = False, box = False)
imagelabel = label(pos = ipos, text='image', yoffset = -2, height = 10, color = color.white, line = False, box = False)
slmeetlabel = label(pos= slmeetpos, text= 'sl meet', yoffset= -2, height = 10, color = color.white, line = False, box = False)

dlslabel = label(pos = ((sx+lx)*.5, oy, oz-15), text = 'DLS: '+ str(dls) +' pc', color = color.white, height = 8, line = False, box = False)
dllabel = label(pos = ((ox+lx)*.5, oy, oz-15), text = 'DL: '+ str(dl) +' pc', color = color.white,height = 8, line = False, box = False)
dslabel = label(pos = ((ox+sx)*.5, oy, oz-25), text = 'DS: '+ str(ds) +' pc', color = color.white, height = 8,line = False, box = False)

