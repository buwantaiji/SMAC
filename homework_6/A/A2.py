import random, math, pylab

def gauss_cut():
    while True:
        x = random.gauss(0.0, 1.0)
        if abs(x) <= 1.0:
            return x

alpha = 0.5
nsteps = 10000000

# A special "partial freezing + direct sampling" algorithm used to get samples of a 2-dimensional distribution pi(x, y).
# Note: essentially, it is a Markov chain Monte Carlo algorithm, as the sample of a iteration depends on the previous sample.
#       However, at the same time, a direct-sampling (concretely, the rejection-sampling, "舍选")approach is used in the coordinate(x or y) that is changed.
samples_x = []
samples_y = []
x, y = 0.0, 0.0
for step in range(nsteps):
    if step % 2 == 0:
        while True:
            x = random.uniform(-1.0, 1.0)
            p = math.exp(-0.5 * x ** 2 - alpha * x ** 4 )
            if random.uniform(0.0, 1.0) < p:
                break
    else:
        while True:
            y = random.uniform(-1.0, 1.0)
            p = math.exp(-0.5 * y ** 2 - alpha * y ** 4 )
            if random.uniform(0.0, 1.0) < p:
                break
    samples_x.append(x)
    samples_y.append(y)
pylab.hexbin(samples_x, samples_y, gridsize=100, bins=1000)
pylab.axis([-1.0, 1.0, -1.0, 1.0])
cb = pylab.colorbar()
pylab.xlabel('x')
pylab.ylabel('y')
pylab.title('A2_1')
#pylab.savefig('plot_A2_1.png')
#pylab.savefig('plot_A2_1.eps')
pylab.show()

# Another version of the algorithm above: propose----one part of pi(x, y); acceptance----the other part of pi(x, y).
samples_x = []
samples_y = []
x, y = 0.0, 0.0
for step in range(nsteps):
    if(step % 2 == 0):
        while(True):
            x = gauss_cut()
            p = math.exp( - alpha * x ** 4 )
            if(random.uniform(0.0, 1.0) < p):
                break
    else:
        while(True):
            y = gauss_cut()
            p = math.exp( - alpha * y ** 4 )
            if(random.uniform(0.0, 1.0) < p):
                break
    samples_x.append(x)
    samples_y.append(y)
pylab.hexbin(samples_x, samples_y, gridsize=100, bins=1000)
pylab.axis([-1.0, 1.0, -1.0, 1.0])
cb = pylab.colorbar()
pylab.xlabel('x')
pylab.ylabel('y')
pylab.title('A2_2')
#pylab.savefig('plot_A2_2.png')
#pylab.savefig('plot_A2_2.eps')
pylab.show()
