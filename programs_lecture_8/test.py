import math, pylab
T = 3.0
mgrid = [i * 0.01 for i in range(-99, 100)]
hgrid = [-2.0 * m ** 2 + 0.5 * T * ((1.0 + m) * math.log(1.0 + m) + (1.0 - m) * math.log(1.0 - m)) for m in mgrid]
pylab.plot(mgrid, hgrid)
pylab.show()
