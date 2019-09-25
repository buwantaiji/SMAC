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

# The algorithm below samples x[0] from the exact probability distribution, which is proportional to rho(x, x, beta).
# Then, it generates the N-1 remaining points of the path using levy construction.
beta = 20.0
N = 80
dtau = beta / N
n_steps = 500000
data = []
for step in range(n_steps):
    xstart = random.gauss( 0.0, 1.0 / math.sqrt(2.0 * math.tanh(0.5 * beta)) )
    x = levy_harmonic_path(xstart, xstart, dtau, N)
#    k = random.randint(0, N - 1)
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
pylab.savefig('plot_B3_beta%s.png' % beta)
pylab.show()
