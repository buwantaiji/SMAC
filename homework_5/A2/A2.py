import math, random, pylab

# This function returns the n-th excited (squared) wave function of 1-dimensional harmonic oscillator: |psi_n(x)|^2
# n = 0 corresponds to the ground state.
def psi_n_square(x, n):
    if n == -1:
        return 0.0
    else:
        psi = [math.exp(-x ** 2 / 2.0) / math.pi ** 0.25]
        psi.append(math.sqrt(2.0) * x * psi[0])
        for k in range(2, n + 1):
            psi.append(math.sqrt(2.0 / k) * x * psi[k - 1] -
                       math.sqrt((k - 1.0) / k) * psi[k - 2])
        return psi[n] ** 2

# This function returns the exact normalized probability distribution pi(x) of 1-dimensional harmonic oscillator at inverse temperature beta, 
#   which is proportional to the diagonal density matrix rho(x, x, beta).
def analytic_prob(x, beta):
    prob = math.exp( - x ** 2 * math.tanh(0.5 * beta)) * math.sqrt( math.tanh(0.5 * beta) / math.pi )
    return prob
# This function returns the analytical probability distribution pi(x) of 1-dimensional harmonic oscillator at inverse temperature beta
#   in the classical limit, which is proportional to exp(- beta * V(x)).
def analytic_prob_classical_lim(x, beta):
    prob = math.exp( - beta * x ** 2 / 2.0) * math.sqrt( beta / (2.0 * math.pi) )
    return prob

beta = 5.0
x = 0.0
n = 0
delta = 0.5
delta_n = 1
iteration_num = 5000000
samples = []
for iteration in range(iteration_num):
    n_new = n + random.choice([-1, 1])
    if random.uniform(0.0, 1.0) < math.exp( - beta * (n_new - n) ) * psi_n_square(x, n_new) / psi_n_square(x, n):
        n = n_new
    x_new = x + random.uniform(-delta, delta)
    if random.uniform(0.0, 1.0) < psi_n_square(x_new, n) / psi_n_square(x, n):
        x = x_new
    samples.append(x)

histo, bin_edges, dummy = pylab.hist(samples, bins=200, normed=True, label="histogram using MCMC sampling method")
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
analytical_result = [analytic_prob(x, beta) for x in bin_centers]
pylab.plot(bin_centers, analytical_result, lw=3, label="analytical value ($\\propto\\rho_{ho}(x, x, \\beta)$)")
analytical_result_classical_lim = [analytic_prob_classical_lim(x, beta) for x in bin_centers]
pylab.plot(bin_centers, analytical_result_classical_lim, lw=3, label="analytical value in classical limit ($\\propto\\mathrm{e}^{-\\beta V(x)}$)")
pylab.title('Probability distribution $\pi(x)$ of a quantum harmonic oscillator at temperature $T=\\frac{1}{\\beta}$\n$\\beta$ = %s' % beta)
pylab.xlabel('$x$', fontsize=18)
pylab.ylabel('$\pi(x)$', fontsize=18)
pylab.legend()
pylab.show()
