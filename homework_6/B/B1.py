import math, random, pylab

def rho_free(x, y, beta):
    return math.exp(-(x - y) ** 2 / (2.0 * beta))

def show_path(x):
    M = len(x)
    path = x + [x[0]]
    y_axis = range(M + 1)
    pylab.plot(path, y_axis, 'bo-')
    pylab.xlim(-5.0, 5.0)
    pylab.xlabel('$x$', fontsize=14)
    pylab.ylabel('$\\tau$', fontsize=14)
    pylab.title("A sampled path with %d slices in PIMC simulation" % M)
    #pylab.savefig("path_snapshot.png")
    pylab.show()

beta = 20.0
N = 80
dtau = beta / N
delta = 1.0
n_steps = 1000000
x = [5.0] * N
data = []
for step in range(n_steps):
    k = random.randint(0, N - 1)
    knext, kprev = (k + 1) % N, (k - 1) % N
    x_new = x[k] + random.uniform(-delta, delta)
    old_weight  = (rho_free(x[knext], x[k], dtau) *
                   rho_free(x[k], x[kprev], dtau) *
                   math.exp(-0.5 * dtau * x[k] ** 2))
    new_weight  = (rho_free(x[knext], x_new, dtau) *
                   rho_free(x_new, x[kprev], dtau) *
                   math.exp(-0.5 * dtau * x_new ** 2))
    if random.uniform(0.0, 1.0) < new_weight / old_weight:
        x[k] = x_new
    if step % N == 0:
#        k = random.randint(0, N - 1)
#        data.append(x[k])
        data += x
show_path(x)

pylab.hist(data, normed=True, bins=400, label='QMC')
list_x = [0.01 * a for a in range (-300, 301)]
list_y = [math.sqrt(math.tanh(beta / 2.0)) / math.sqrt(math.pi) * \
          math.exp(-x ** 2 * math.tanh(beta / 2.0)) for x in list_x]
pylab.plot(list_x, list_y, label='analytic')
pylab.legend()
pylab.xlabel('$x$')
pylab.ylabel('$\\pi(x)$ (normalized)')
pylab.title('naive_harmonic_path (beta=%s, N=%i)' % (beta, N))
pylab.xlim(-3.0, 3.0)
#pylab.savefig('plot_B1_beta%s.png' % beta)
pylab.show()
