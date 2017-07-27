from visual import *
import numpy as np

ox = -80.0
oy = 0.0
oz = -40.0
opos = vector(ox, oy, oz)
#observer.pos = (ox, oy)

lx = 0.0
ly = 0.0
lz = 40.0
lpos = vector(lx, ly, lz)
#lens.pos = (lx, ly)


sx = 80.0
sy = 0.0
sz = 40.0
spos = vector(sx, sy, sz)
#source.pos = (sx, sy)

ds = sx-ox
dl = lx-ox
dls = sx-lx

#rt = paths.rectangle(pos=(0,0), height= 200, width = 200)
#curve(pos=rt.pos)
ols = points(pos=[spos, lpos, opos], size = 10, color = color.red)


osline = paths.line(start=spos, end =opos, np = 20)
curve(pos=osline.pos)
olline = paths.line(start=opos, end=lpos, np=20)
curve(pos=olline.pos)

dsline = paths.line(start=opos, end=(sx, oy, oz), np = 20)
curve(pos=dsline.pos)
dlline = paths.line(start=(ox, oy, oz-4), end=(lx-1, oy, oz-4), np=20)
curve(pos=dlline.pos)
dlsline = paths.line(start=(lx+1, oy, oz-4), end=(sx, oy, oz-4), np=20)
curve(pos=dlsline.pos)

sperpline = paths.line(start=spos, end=(sx, oy, oz), np=20)
curve(pos=sperpline.pos)
lperpline = paths.line(start=lpos, end=(lx, oy, oz), np = 50)
curve(pos=lperpline.pos)

slline = paths.line(start=lpos, end=spos)
curve(pos=slline.pos)

ok = lz-oz

theta = atan2(ok,dl)
psz = (tan(theta))*(ds)
pspos = vector(sx, sy, psz-lz)
ps = points(pos=pspos, size = 10, color = color.cyan)
opsline = paths.line(start=opos, end=pspos, np = 20)
curve(pos=opsline.pos)

idk = psz-oz
t2 = atan2(idk, ds)*180/pi
print(ds)
