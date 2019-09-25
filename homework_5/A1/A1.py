import random, math, pylab

def psi0_sq(x):
    psisq = math.exp( - x ** 2 ) / math.sqrt(math.pi)
    return psisq

x = 0.0
delta = 0.5
samples = []
for k in range(5000000):
    x_new = x + random.uniform(-delta, delta)
    if random.uniform(0.0, 1.0) < psi0_sq(x_new) / psi0_sq(x):
        x = x_new 
    samples.append(x)

histo, bin_edges, dummy = pylab.hist(samples, bins=200, normed=True, label="histogram using MCMC sampling method")
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
analytical_result = [psi0_sq(x) for x in bin_centers]
pylab.plot(bin_centers, analytical_result, "r-", lw=3, label="analytical value(=$|\psi_0(x)|^2$)")
pylab.title('Probability distribution $\pi(x)$ of a quantum harmonic oscillator at temperature T=0')
pylab.xlabel('$x$', fontsize=18)
pylab.ylabel('$\pi(x)$', fontsize=18)
pylab.legend()
pylab.savefig('pi_x_T_0.png')
pylab.show()
