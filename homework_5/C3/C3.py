import math, numpy, pylab

# Free off-diagonal density matrix
def rho_free(x, xp, beta):
    return (math.exp(-(x - xp) ** 2 / (2.0 * beta)) /
            math.sqrt(2.0 * math.pi * beta))

# The anharmonic potential, including cubic and quartic terms.
def anharmonic(x, cubic, quartic):
    V = x ** 2 / 2.0 + cubic * x ** 3 + quartic * x ** 4
    return V

# Anharmonic density matrix in the Trotter approximation (returns the full matrix)
def rho_anharmonic_trotter(grid, beta, cubic, quartic):
    return numpy.array([[rho_free(x, xp, beta) * \
                         numpy.exp( -0.5 * beta * (anharmonic(x, cubic, quartic) + anharmonic(xp, cubic, quartic)) ) \
                         for x in grid] for xp in grid])

def Energy_pert(n, cubic, quartic):
    return n + 0.5 - 15.0 / 4.0 * cubic **2 * (n ** 2 + n + 11.0 / 30.0) \
         + 3.0 / 2.0 * quartic * (n ** 2 + n + 1.0 / 2.0)

def Z_pert(cubic, quartic, beta, n_max):
    Z = sum(math.exp(-beta * Energy_pert(n, cubic, quartic)) for n in range(n_max + 1))
    return Z

quartics = [0.0, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5]

x_max = 5.0
nx = 1000
dx = 2.0 * x_max / (nx - 1)
x = [i * dx for i in range(-(nx - 1) / 2, nx / 2 + 1)]

# Plot the harmonic and anharmonic potential here.
for i in range(len(quartics)):
    quartic = quartics[i]
    cubic = - quartic
    potential= [anharmonic(x[j], cubic, quartic) for j in range(nx + 1)]
    pylab.plot(x, potential, lw=2, label="quartic = %s" % quartic)
pylab.title("Various anharmonic potential\n"
                + "$V(x) = \\frac{1}{2} x^2 + \mathrm{cubic} x^3 + \mathrm{quartic} x^4.$\n"
                + "cubic = - quartic")
pylab.xlabel('$x$', fontsize=18)
pylab.ylabel('$V(x)$', fontsize=18)
pylab.legend()
pylab.show()

# Compute single particle probability distribution pi(x) using matrix squaring.
# In addition, in the regime of small quartic and cubic, compute the computed partition function with the analytical perturbative result.
beta_init = 2.0 ** (-6)                                       # initial value of beta (power of 2)
beta      = 2.0 ** 1                                          # actual value of beta (power of 2)
n_max = 20
for i in range(len(quartics)):
    quartic = quartics[i]
    cubic = - quartic
    rho= rho_anharmonic_trotter(x, beta_init, cubic, quartic)    # density matrix at initial beta
    beta_iter = beta_init
    while beta_iter < beta:
        rho= numpy.dot(rho, rho)
        rho*= dx
        beta_iter *= 2.0
        #print 'beta: %s -> %s' % (beta_iter / 2.0, beta_iter)
    Z= sum(rho[j, j] for j in range(nx + 1)) * dx
    Z_perturbation = Z_pert(cubic, quartic, beta, n_max)
    print "quartic = %05s:\t\tZ = %s\t\tZ_perturbation = %s\n" % (quartic, Z, Z_perturbation)
    pi_of_x= [rho[j, j] / Z for j in range(nx + 1)]
    pylab.plot(x, pi_of_x, lw=2, label="quartic = %s" % quartic)
pylab.xlim(-2.0, 2.0)
pylab.title("Probability distribution $\pi(x)$ of a particle in various anharmonic potential at temperature $T=\\frac{1}{\\beta}$\n"
                + "$V(x) = \\frac{1}{2} x^2 + \mathrm{cubic} x^3 + \mathrm{quartic} x^4.$\n"
                + "cubic = - quartic\n"
                + "$\\beta$ = %s" % beta)
pylab.xlabel('$x$', fontsize=18)
pylab.ylabel('$\pi(x)$', fontsize=18)
pylab.legend()
pylab.show()
