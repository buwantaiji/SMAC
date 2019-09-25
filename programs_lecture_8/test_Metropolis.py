import math, random, pylab

L = 128
N = L * L
nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
            (i // L) * L + (i - 1) % L, (i - L) % N)
                                    for i in range(N)}
B = -0.01
T = 1.0
beta = 1.0 / T
S = [random.choice([1, -1]) for site in range(N)]
m = sum(S) / float(N)
m_samples = []

nsteps = N * 10000
nsteps_per_record = N * 100
for step in range(nsteps):
    site = random.randint(0, N - 1)
    delta_E = 2.0 * S[site] * ( sum(S[nn] for nn in nbr[site]) + B )
    if(random.uniform(0.0, 1.0) < math.exp(- beta * delta_E)):
        m -= 2 * S[site] / float(N)
        S[site] *= -1
    if(step % nsteps_per_record == 0):
        print "%d: m = %s" % (step / nsteps_per_record, m)
    m_samples.append(m)

pylab.plot(range(nsteps), m_samples)
pylab.ylim(-1.0, 1.0)
pylab.show()
