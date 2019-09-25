import random, math, pylab

def energy(S, N, nbr):
    E = 0.0
    for k in range(N):
        E -=  S[k] * sum(S[nn] for nn in nbr[k])
    return 0.5 * E / N

L = 32
N = L * L
nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
            (i // L) * L + (i - 1) % L, (i - L) % N)
                                    for i in range(N)}

T = 2.27
p  = 1.0 - math.exp(-2.0 / T)
nsteps = 10000
S = [random.choice([1, -1]) for k in range(N)]
E_records = [energy(S, N, nbr)]
M_records = [sum(S) / float(N)]
for step in range(nsteps):
    k = random.randint(0, N - 1)
    Pocket, Cluster = [k], [k]
    while Pocket != []:
        j = random.choice(Pocket)
        for l in nbr[j]:
            if S[l] == S[j] and l not in Cluster \
                   and random.uniform(0.0, 1.0) < p:
                Pocket.append(l)
                Cluster.append(l)
        Pocket.remove(j)
    for j in Cluster:
        S[j] *= -1
    E_records.append(energy(S, N, nbr))
    M_records.append(sum(S) / float(N))
    if(step % 100 == 0):
        print "%d: M = %s" % (step / 100, M_records[-1])
print 'mean energy per spin:', sum(E_records) / float(len(E_records))

pylab.plot(range(nsteps + 1), E_records, "g")
pylab.xlabel("steps")
pylab.ylabel("$\\frac{E}{N}$")
pylab.title("Energy per spin of 2-dimentional Ising model at temperature $T$\n"
            + "Square lattice: $%d\\times%d$, periodic boundary condition\n" % (L, L)
            + "$T=%s$" % T)
pylab.show()
pylab.plot(range(nsteps + 1), M_records, "ro")
pylab.xlabel("steps")
pylab.ylabel("$\\frac{M}{N}$")
pylab.ylim(- 1.0, 1.0)
pylab.title("Magnetization per spin of 2-dimentional Ising model at temperature $T$\n"
            + "Square lattice: $%d\\times%d$, periodic boundary condition\n" % (L, L)
            + "$T=%s$" % T)
pylab.show()
