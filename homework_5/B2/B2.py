import math, random, pylab

def read_file(filename):
    list_x = []
    list_y = []
    with open(filename) as f:
        for line in f:
            x, y = line.split()
            list_x.append(float(x))
            list_y.append(float(y))
    f.close()
    return list_x, list_y

def rho_free(x, y, beta):    # free off-diagonal density matrix
    return math.exp(-(x - y) ** 2 / (2.0 * beta)) 

beta = 4.0
N = 10                                            # number of slices
dtau = beta / N
delta = 1.0                                       # maximum displacement on one slice
n_steps = 10000000                                 # number of Monte Carlo steps
#n_steps = 100
steps_per_sample = 10
samples = []
x = [0.0] * N                                     # initial path
for step in range(n_steps):
    k = random.randint(0, N - 1)                  # random slice
    knext, kprev = (k + 1) % N, (k - 1) % N       # next/previous slices
    x_new = x[k] + random.uniform(-delta, delta)  # new position at slice k
    old_weight  = (rho_free(x[knext], x[k], dtau) *
                   rho_free(x[k], x[kprev], dtau) *
                   math.exp(-0.5 * dtau * x[k] ** 2))
    new_weight  = (rho_free(x[knext], x_new, dtau) *
                   rho_free(x_new, x[kprev], dtau) *
                   math.exp(-0.5 * dtau * x_new ** 2))
    if random.uniform(0.0, 1.0) < new_weight / old_weight:
        x[k] = x_new
#    print "iter %3d:" % (step + 1), x
    if(step % steps_per_sample == 0):
        samples += x

pylab.hist(samples, bins=200, normed=True, label="histogram using PIMC sampling method")

matrixsquaring_data_filename = "pi_x_matrixsquaring_beta_4.0.dat"
position, pi_x = read_file(matrixsquaring_data_filename)
pylab.plot(position, pi_x, lw=2, label="direct matrix squaring of density matrices")

pylab.title("Probability distribution $\pi(x)$(i.e. the normalized 'single particle density') of a quantum harmonic oscillator\n$\\beta$ = %s" % beta)
pylab.xlabel('$x$', fontsize=18)
pylab.ylabel('$\pi(x)$', fontsize=18)
pylab.legend()
pylab.show()
