import math, numpy, pylab

# This function returns the exact normalized probability distribution pi(x) of 1-dimensional harmonic oscillator at inverse temperature beta, 
#   which is proportional to the diagonal density matrix rho(x, x, beta).
def analytic_prob(x, beta):
    prob = math.exp( - x ** 2 * math.tanh(0.5 * beta)) * math.sqrt( math.tanh(0.5 * beta) / math.pi )
    return prob

# Free off-diagonal density matrix
def rho_free(x, xp, beta):
    return (math.exp(-(x - xp) ** 2 / (2.0 * beta)) /
            math.sqrt(2.0 * math.pi * beta))

# Harmonic density matrix in the Trotter approximation (returns the full matrix)
def rho_harmonic_trotter(grid, beta):
    return numpy.array([[rho_free(x, xp, beta) * \
                         numpy.exp(-0.5 * beta * 0.5 * (x ** 2 + xp ** 2)) \
                         for x in grid] for xp in grid])

x_max = 5.0
nx = 1000
dx = 2.0 * x_max / (nx - 1)
x = [i * dx for i in range(-(nx - 1) / 2, nx / 2 + 1)]
beta_tmp = 2.0 ** (-6)                   # initial value of beta (power of 2)
beta     = 2.0 ** 2                      # actual value of beta (power of 2)
rho = rho_harmonic_trotter(x, beta_tmp)  # density matrix at initial beta
while beta_tmp < beta:
    rho = numpy.dot(rho, rho)
    rho *= dx
    beta_tmp *= 2.0
    print 'beta: %s -> %s' % (beta_tmp / 2.0, beta_tmp)

Z = sum(rho[j, j] for j in range(nx + 1)) * dx
pi_of_x = [rho[j, j] / Z for j in range(nx + 1)]
#f = open('pi_x_matrixsquaring_beta_' + str(beta) + '.dat', 'w')
#for j in range(nx + 1):
    #f.write(str(x[j]) + ' ' + str(rho[j, j] / Z) + '\n')
#f.close()

analytical_result = [analytic_prob(pos, beta) for pos in x]
pylab.plot(x, analytical_result, lw=2, label="analytical value ($\\propto\\rho_{ho}(x, x, \\beta)$)")
pylab.plot(x, pi_of_x, lw=2, label="matrix squaring method")
#pylab.xlim(-x_max, x_max)
pylab.title('Probability distribution $\pi(x)$ of a quantum harmonic oscillator at temperature $T=\\frac{1}{\\beta}$\n$\\beta$ = %s' % beta)
pylab.xlabel('$x$', fontsize=18)
pylab.ylabel('$\pi(x)$', fontsize=18)
pylab.legend()
pylab.show()
