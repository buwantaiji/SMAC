import random, math, pylab

def energy(S, N, nbr):
    E = 0.0
    for k in range(N):
        E -=  S[k] * sum(S[nn] for nn in nbr[k])
    return 0.5 * E / N

L = 32
N = L * L
nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
            (i // L) * L + (i - 1) % L, (i - L) % N) \
                                    for i in range(N)}

T = 2.2
S = [random.choice([1, -1]) for k in range(N)]
nsteps = N * 100000
nsteps_per_record = N * 100
beta = 1.0 / T
Energy = energy(S, N, nbr)
Magnetization = sum(S) / float(N)
E_records = []
M_records = []
for step in xrange(nsteps):
    k = random.randint(0, N - 1)
    delta_E = 2.0 * S[k] * sum(S[nn] for nn in nbr[k])
    delta_M = - 2 * S[k]
    if random.uniform(0.0, 1.0) < math.exp(-beta * delta_E):
        S[k] *= -1
        Energy += delta_E / float(N)
        Magnetization += delta_M / float(N)
    if(step % (nsteps_per_record / 100) == 0):
        E_records.append(Energy)
        M_records.append(Magnetization)
    if(step % nsteps_per_record == 0):
        print "%d: M = %s" % (step / nsteps_per_record, Magnetization)
   
print 'mean energy per spin:', sum(E_records) / float(len(E_records))

pylab.plot([step / float(100) for step in range(1, len(E_records) + 1)], E_records, "g")
pylab.xlabel("steps\n(in unit of $N * 100$)")
pylab.ylabel("$\\frac{E}{N}$")
pylab.title("Energy per spin of 2-dimentional Ising model at temperature $T$\n"
            + "Square lattice: $%d\\times%d$, periodic boundary condition\n" % (L, L)
            + "$T=%s$" % T)
pylab.show()
pylab.plot([step / float(100) for step in range(1, len(M_records) + 1)], M_records, "r")
pylab.xlabel("steps\n(in unit of $N * 100$)")
pylab.ylabel("$\\frac{M}{N}$")
pylab.ylim(- 1.0, 1.0)
pylab.title("Magnetization per spin of 2-dimentional Ising model at temperature $T$\n"
            + "Square lattice: $%d\\times%d$, periodic boundary condition\n" % (L, L)
            + "$T=%s$" % T)
pylab.show()
