import random, pylab

N = 24
L = 30.0
sigma = 0.5
n_runs = 80000
data = []
for run in range(n_runs):
    Lprime = L - 2.0 * sigma
    y_sorted = [random.uniform(0, Lprime - 2.0 * (N - 1) * sigma) for k in range(N - 1)]
    y_sorted.sort()
    sample = [y_sorted[k] + (2.0 * k + 1.0) * sigma for k in range(N - 1)] + [L - sigma]
    # Take care of the periodic boundary condition.
    # When sampling, one of the N particles is "free" to be any position between 0 and L. In fact:
    # The configuration integral(spatial part of partition function) for two kinds of boundary conditions:
    #   "Rigid wall" boundary condition: (L - 2 * N * sigma)^N
    #   Periodic boundary condition: L(L - 2 * N * sigma)^(N - 1)
    shift = random.uniform(0, L)
    data += [(y + shift) % L for y in sample]
pylab.title('Density of %i clothes-pins ($\sigma$=%s) on a line of length L=%s' % (N, sigma, L))
pylab.xlabel('$x$', fontsize=14)
pylab.ylabel('$\pi(x)$', fontsize=14)
pylab.title('Density profile $\pi(x)$ for N=%i, $\sigma$=%.2f, L=%.1f' % (N, sigma, L))
pylab.hist(data, bins=200, normed=True)

# Analytical approach for (normalized) singleparticle density under periodic boundary condition.
# Notice that under periodic boundary condition, the system is homogeneous, and the result we obtained is exactly uniform!
xr = pylab.linspace(0.0, L, 201)
yr = [1.0 / L for x in xr]
pylab.plot(xr, yr, 'red', linewidth=2.0)

pylab.savefig('plot-pins_noreject_periodic-N%04i-L%.1f-density.png' % (N, L))
pylab.show()
