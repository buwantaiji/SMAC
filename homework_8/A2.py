import os, random, math, pylab

def energy(S, N, nbr):
    E = 0.0
    for k in range(N):
        E -=  S[k] * sum(S[nn] for nn in nbr[k])
    return 0.5 * E

L = 128
N = L * L
nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
            (i // L) * L + (i - 1) % L, (i - L) % N) \
                                    for i in range(N)}
T = 2.27
beta = 1.0 / T

filename = "%d*%d/snapshots/configuration_profile_L_%d_T_%s.txt" % (L, L, L, T)
if os.path.isfile(filename):
    f = open(filename, 'r')
    S = []
    for line in f:
        S.append(int(line))
    f.close()
    print 'Starting from file', filename
else:
    S = [random.choice([1, -1]) for k in range(N)]
    print 'Starting from a random configuration'

nsteps = N * 50000
nsteps_per_record = N * 100
Energy = energy(S, N, nbr)
Magnetization = sum(S) / float(N)
for step in xrange(nsteps):
    k = random.randint(0, N - 1)
    delta_E = 2.0 * S[k] * sum(S[nn] for nn in nbr[k])
    delta_M = - 2 * S[k]
    if random.uniform(0.0, 1.0) < math.exp(-beta * delta_E):
        S[k] *= -1
        Energy += delta_E
        Magnetization += delta_M / float(N)
    if(step % nsteps_per_record == 0):
        print "%d: M = %s" % (step / nsteps_per_record, Magnetization)
print "The final configuration: M = %s" % Magnetization

f = open(filename, 'w')
for a in S:
   f.write(str(a) + '\n')
f.close()

# conf[x][y] = S[k], where k = y * L + x.
conf = [[S[y * L + x] for y in range(L)] for x in range(L)]

pylab.imshow(conf, extent=[0, L, 0, L], interpolation='nearest')
pylab.set_cmap('hot')
pylab.title("Typical equilibrium configuration of 2-dimensional Ising model at temperature T\n"
            + "Square lattice: $%d\\times%d$, periodic boundary condition\n" % (L, L)
            + "$T=%s$" % T)
#pylab.savefig("%d*%d/snapshots/L_%d_T_%s_snapshot.png" % (L, L, L, T))
pylab.show()
