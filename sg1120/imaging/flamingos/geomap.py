#!/usr/bin/python
from numarray import *
import math, scipy
from scipy.optimize import leastsq
import Numeric

infile = sys.argv[1]
order = int(sys.argv[2])
outfile = sys.argv[3]
center = float(sys.argv[4])
terms = 0
for j in range(order+2):
   terms+=j
f = open(infile,'rb')
a = f.read()
f.close()
a = a.strip()
a = a.split('\n')
xin = []
xout = []
yin = []
yout = []
for j in a:
   temp=j.split()
   xout.append(float(temp[0]))
   yout.append(float(temp[1]))
   xin.append(float(temp[2]))
   yin.append(float(temp[3]))
xout = array(xout)-center
yout = array(yout)-center
xin = array(xin)-center
yin = array(yin)-center

def residuals(p, x, y, out, order):
   f = zeros(len(x))+0.
   n = 0
   for j in range(order+1):
      for l in range(j+1):
        f+=p[n]*x**(j-l)*y**l
        n+=1
   err = out - f
   err = Numeric.array(err.tolist())
   return err

p0x = zeros(terms)
p0y = zeros(terms)
p0x[1] = 1
p0y[2] = 1
p0x = Numeric.array(p0x.tolist())
p0y = Numeric.array(p0y.tolist())

xlsq = leastsq(residuals, p0x, args=(xin, yin, xout, order))
ylsq = leastsq(residuals, p0y, args=(xin, yin, yout, order))

f = open(outfile,'wb')
f.write('poly ' + str(order) + '\n')
for j in xlsq[0]:
   f.write(str(j)+'\n')
f.write('\n')
for j in ylsq[0]:
   f.write(str(j)+'\n')
f.close()
