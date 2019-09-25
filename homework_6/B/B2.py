import math, random, pylab
def levy_harmonic_path(xstart, xend, dtau, N):
    x = [xstart]
    for k in range(1, N):
        dtau_prime = (N - k) * dtau
        Ups1 = 1.0 / math.tanh(dtau) + \
               1.0 / math.tanh(dtau_prime)
        Ups2 = x[k - 1] / math.sinh(dtau) + \
               xend / math.sinh(dtau_prime)
        x.append(random.gauss(Ups2 / Ups1, \
               1.0 / math.sqrt(Ups1)))
    return x

# The algorithm below generates the intermidiate N-1 points of a path using levy construction.
# However, the "starting point" x[0] is arbitrarily initialized, where the "partial freezing" part enters.
beta = 20.0
N = 2
dtau = beta / N
n_steps = 5000000
data = []
xstart = 5.0
for step in range(n_steps):
    x = levy_harmonic_path(xstart, xstart, dtau, N)
#    N_cut = N / 2
#    x = x[N_cut:] + x[:N_cut]
    k = random.randint(0, N - 1)
    xstart = x[k]
#    data.append(x[k])
    data += x

pylab.hist(data, normed=True, bins=400, label='QMC')
list_x = [0.1 * a for a in range (-30, 31)]
list_y = [math.sqrt(math.tanh(beta / 2.0)) / math.sqrt(math.pi) * \
          math.exp(-x ** 2 * math.tanh(beta / 2.0)) for x in list_x]
pylab.plot(list_x, list_y, label='analytic')
pylab.legend()
pylab.xlabel('$x$')
pylab.ylabel('$\\pi(x)$ (normalized)')
pylab.title('naive_harmonic_path (beta=%s, N=%i)' % (beta, N))
pylab.xlim(-2, 2)
pylab.savefig('plot_B2_beta%s.png' % beta)
pylab.show()
