import random, math, pylab

def gauss_cut():
    while True:
        x = random.gauss(0.0, 1.0)
        if abs(x) <= 1.0:
            return x

alpha = 1.0
nsteps = 50000000

# The original Metropolis MC algorithm: propose----uniform; acceptance----full pi(x, y).
# This is the generic Markov Chain Monte Carlo sampling method, without use of rejection-sampling or any other kinds of direct sampling methods.
samples_x = []
samples_y = []
x, y = 0.0, 0.0
for step in range(nsteps):
    xnew = random.uniform(-1.0, 1.0)
    ynew = random.uniform(-1.0, 1.0)
    exp_new = - 0.5 * (xnew ** 2 + ynew ** 2) - alpha * (xnew ** 4 + ynew ** 4)
    exp_old = - 0.5 * (x ** 2 + y ** 2) - alpha * (x ** 4 + y ** 4)
    if random.uniform(0.0, 1.0) < math.exp(exp_new - exp_old):
        x = xnew
        y = ynew
    samples_x.append(x)
    samples_y.append(y)
pylab.hexbin(samples_x, samples_y, gridsize=100, bins=1000)
pylab.axis([-1.0, 1.0, -1.0, 1.0])
cb = pylab.colorbar()
pylab.xlabel('x')
pylab.ylabel('y')
pylab.title('A3_1')
pylab.savefig('plot_A3_1_alpha_%s.png' % alpha)
pylab.savefig('plot_A3_1_alpha_%s.eps' % alpha)
pylab.show()

# A modified version of the Metropolis MC algorithm: propose----one part of pi(x, y); acceptance----the other part of pi(x, y).
samples_x = []
samples_y = []
x, y = 0.0, 0.0
for step in range(nsteps):
    x_new, y_new = gauss_cut(), gauss_cut()
    exp_new = - alpha * (x_new ** 4 + y_new ** 4)
    exp_old = - alpha * (x ** 4 + y ** 4)
    if random.uniform(0.0, 1.0) < math.exp(exp_new - exp_old):
        x, y = x_new, y_new
    samples_x.append(x)
    samples_y.append(y)
pylab.hexbin(samples_x, samples_y, gridsize=100, bins=1000)
pylab.axis([-1.0, 1.0, -1.0, 1.0])
cb = pylab.colorbar()
pylab.xlabel('x')
pylab.ylabel('y')
pylab.title('A3_2')
pylab.savefig('plot_A3_2_alpha_%s.png' % alpha)
pylab.savefig('plot_A3_2_alpha_%s.eps' % alpha)
pylab.show()
