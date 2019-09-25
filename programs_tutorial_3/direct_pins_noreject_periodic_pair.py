# This program computes the (normalized) pair correlation function for 1-dimensional "hard sphere" system.
# Both direct sampling MC simulation(with no rejections) and analytical results are obtained.
import random, pylab

def binomialCoeff(n, k):
    result = 1
    for i in range(1, k+1):
        result = result * (n-i+1) / i
    return result

def Z(N, L, sigma):
    freespace = L - 2.0 * N * sigma
    if freespace > 0.0:
        result = freespace ** N
    else:
        result = 0.0
    return result

def pi(x, N, L, sigma):
    tot = 0.
    for k in range(0, N):
        Z1 = Z(k, x - sigma, sigma)
        Z2 = Z(N - k - 1, L - x - sigma, sigma)
        tot += binomialCoeff( N - 1, k) * Z1 * Z2
    Ztotal = Z(N, L, sigma)
    return tot / Ztotal

def dist(x1, x2, L):
    d_x = abs(x1 - x2) 
    return min(d_x, L - d_x)

#N = 450
#L = 500.0
#sigma = 0.5
#density = N * 2.0 * sigma / L
N = 75
L = 200.0
density = 0.92
sigma = density * L / (2.0 * N)
n_runs = 8000
#x_max = 30.0  # maximum of the histogram range
x_max = L / 2.0

# The direct sampling MC approach, which can be non-rejection for 1-dimensional case.
data, pair_corr = [], []
for run in range(n_runs):
    Lprime = L - 2.0 * sigma
    y_sorted = [random.uniform(0, Lprime - 2.0 * (N - 1.0) * sigma) for k in xrange(N - 1)]
    y_sorted.sort()
    sample = [y_sorted[k] + (2.0 * k + 1.0) * sigma for k in xrange(N - 1) ] + [L - sigma]
    # Think carefully about the sampling here, which is exactly sampling the (normalized) probability distribution pi(x, x') = pi(|x - x'|)
    # This distribution pi(x) is related to the usual pair correlation function g(x) by: g(x) = L * pi(x). (Generally, in three dimension, we have g(r) = V * pi(r).)
    pair_corr += [dist(sample[i], sample[j], L) for i in xrange(N) for j in xrange(i)]
histo, bins, patches = pylab.hist(pair_corr, bins=800, normed=True)

# The analytical approach. The figure we finally obtain will indicate that the MC simulation and analytical results match very well!!!
distance = pylab.linspace(sigma, x_max, 801)
xr = [L - sigma - d for d in distance]
# Note: the factor "2.0" is included to account for the fact that the distance d only take one half of the whole length L due to periodic boundary condition.
#   Thus, the factor "2.0" is included to fix the probability normalization.
pi_analytic = [2.0 * pi(x, N - 1, L - 2.0 * sigma, sigma) for x in xr]
pi_no_correlation = [1.0 / x_max for d in distance]
pylab.plot(distance, pi_analytic, 'red', linewidth=2.0)
pylab.plot(distance, pi_no_correlation, 'black', linewidth=1.0)

pylab.xlim(0.0, x_max)
pylab.title('(Normalized) Pair-correlation function $\pi(x,y)$\nN=%i, $\sigma$=%.2f, L=%.1f, density=%.2f' % (N, sigma, L, density))
pylab.xlabel('$|x-y|$', fontsize=14)
pylab.ylabel('$\pi(|x-y|)$', fontsize=14)
pylab.savefig('pair_correlation/pair_correlation_function(normalized)_eta_%f.eps' % density)
pylab.show()
pylab.clf()
asymptotic_val = 1.0 / (L / 2.0)   # asymptotic (uncorrelated) value of the pair correlation function
#pylab.semilogy(bins[:-1], [abs(y - asymptotic_val) for y in histo])
pylab.plot(distance, [pi - asymptotic_val for pi in pi_analytic])
pylab.xlim(0.0, x_max)
pylab.title('Deviation of $\pi(x,y)$ from its asymptotic (uncorrelated) value\nN=%i, $\sigma$=%.2f, L=%.1f, density=%.2f' % (N, sigma, L, density))
pylab.xlabel('$|x-y|$', fontsize=14)
pylab.ylabel('$\pi(|x-y|)-\pi_\mathrm{asympt}$\n$\propto$\n$cov(n(x), n(y)) = <n(x)n(y)> - <n(x)><n(y)>$', fontsize=14)
#pylab.ylabel('$|\pi(|x-y|)-\pi_\mathrm{asympt}|$', fontsize=14)
#pylab.savefig('pair_correlation/correlation_deviation_eta_%f.eps' % density)
pylab.show()
