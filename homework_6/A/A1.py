import random, math, pylab

def gauss_cut():
    while True:
        x = random.gauss(0.0, 1.0)
        if abs(x) <= 1.0:
            return x

alpha = 0.5
nsamples = 10000000

# This is a pure rejection-sampling algorithm(""), which is used to get correct samples of a 2-dimensional probability distribution.
samples_x = []
samples_y = []
for sample in xrange(nsamples):
    while True:
        x = random.uniform(-1.0, 1.0)
        y = random.uniform(-1.0, 1.0)
        p = math.exp(-0.5 * (x ** 2 + y ** 2) - alpha * (x ** 4 + y ** 4))
        if random.uniform(0.0, 1.0) < p:
            break
    samples_x.append(x)
    samples_y.append(y)
pylab.hexbin(samples_x, samples_y, gridsize=100, bins=1000)
pylab.axis([-1.0, 1.0, -1.0, 1.0])
cb = pylab.colorbar()
pylab.xlabel('x')
pylab.ylabel('y')
pylab.title('A1_1')
#pylab.savefig('plot_A1_1.png')
pylab.show()

# Another version of the rejection-sampling algorithm above: propose: one part of pi(x, y); acceptance: the other part of pi(x, y).
samples_x = []
samples_y = []
for sample in xrange(nsamples):
    while True:
        x = gauss_cut()
        y = gauss_cut()
        p = math.exp( - alpha * (x ** 4 + y ** 4) )
        if random.uniform(0.0, 1.0) < p:
            break
    samples_x.append(x)
    samples_y.append(y)
pylab.hexbin(samples_x, samples_y, gridsize=100, bins=1000)
pylab.axis([-1.0, 1.0, -1.0, 1.0])
cb = pylab.colorbar()
pylab.xlabel('x')
pylab.ylabel('y')
pylab.title('A1_2')
#pylab.savefig('plot_A1_2.png')
pylab.show()
