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

# The anharmonic potential, including cubic and quartic terms.
def anharmonic(x, cubic, quartic):
    V = x ** 2 / 2.0 + cubic * x ** 3 + quartic * x ** 4
    return V

# Harmonic density matrix in the Trotter approximation (returns the full matrix)
def rho_harmonic_trotter(grid, beta):
    return numpy.array([[rho_free(x, xp, beta) * \
                         numpy.exp(-0.5 * beta * 0.5 * (x ** 2 + xp ** 2)) \
                         for x in grid] for xp in grid])
# Anharmonic density matrix in the Trotter approximation (returns the full matrix)
def rho_anharmonic_trotter(grid, beta, cubic, quartic):
    return numpy.array([[rho_free(x, xp, beta) * \
                         numpy.exp( -0.5 * beta * (anharmonic(x, cubic, quartic) + anharmonic(xp, cubic, quartic)) ) \
                         for x in grid] for xp in grid])

quartic = 1.0
cubic = - quartic

x_max = 5.0
nx = 1000
dx = 2.0 * x_max / (nx - 1)
x = [i * dx for i in range(-(nx - 1) / 2, nx / 2 + 1)]

# Plot the harmonic and anharmonic potential here.
potential_harmonic   = [0.5 * x[j] ** 2 for j in range(nx + 1)]
potential_anharmonic = [anharmonic(x[j], cubic, quartic) for j in range(nx + 1)]
pylab.plot(x, potential_harmonic,   lw=2, label="harmonic")
pylab.plot(x, potential_anharmonic, lw=2, label="anharmonic(cubic = %s, quartic = %s)" % (cubic, quartic))
pylab.title("Harmonic and anharmonic potential")
pylab.xlabel('$x$', fontsize=18)
pylab.ylabel('$V(x)$', fontsize=18)
pylab.legend()
pylab.show()

beta_tmp = 2.0 ** (-6)                                       # initial value of beta (power of 2)
beta     = 2.0 ** 2                                          # actual value of beta (power of 2)
rho_harmonic = rho_harmonic_trotter(x, beta_tmp)
rho_anharmonic = rho_anharmonic_trotter(x, beta_tmp, cubic, quartic)    # density matrix at initial beta
while beta_tmp < beta:
    rho_harmonic = numpy.dot(rho_harmonic, rho_harmonic)
    rho_anharmonic = numpy.dot(rho_anharmonic, rho_anharmonic)
    rho_harmonic *= dx
    rho_anharmonic *= dx
    beta_tmp *= 2.0
    print 'beta: %s -> %s' % (beta_tmp / 2.0, beta_tmp)

Z_harmonic   = sum(rho_harmonic[j, j] for j in range(nx + 1)) * dx
Z_anharmonic = sum(rho_anharmonic[j, j] for j in range(nx + 1)) * dx
pi_of_x_harmonic   = [rho_harmonic[j, j] / Z_harmonic for j in range(nx + 1)]
pi_of_x_anharmonic = [rho_anharmonic[j, j] / Z_anharmonic for j in range(nx + 1)]
#f = open('pi_x_anharmonic_matrixsquaring_beta_' + str(beta) + '.dat', 'w')
#for j in range(nx + 1):
    #f.write(str(x[j]) + ' ' + str(rho_anharmonic[j, j] / Z_anharmonic) + '\n')
#f.close()

#analytical_result = [analytic_prob(pos, beta) for pos in x]
#pylab.plot(x, analytical_result, lw=2, label="analytical value ($\\propto\\rho_{ho}(x, x, \\beta)$)")
pylab.plot(x, pi_of_x_harmonic,   lw=2, label="harmonic")
pylab.plot(x, pi_of_x_anharmonic, lw=2, label="anharmonic(cubic = %s, quartic = %s)" % (cubic, quartic))
pylab.xlim(-2.0, 2.0)
pylab.title('Probability distribution $\pi(x)$ of a particle in a harmonic/anharmonic potential at temperature $T=\\frac{1}{\\beta}$\n$\\beta$ = %s' % beta)
pylab.xlabel('$x$', fontsize=18)
pylab.ylabel('$\pi(x)$', fontsize=18)
pylab.legend()
pylab.show()
