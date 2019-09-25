import os, random, math

def energy(S, N, nbr):
    E = 0.0
    for k in range(N):
        E -=  S[k] * sum(S[nn] for nn in nbr[k])
    return 0.5 * E

L = 32
N = L * L
nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
            (i // L) * L + (i - 1) % L, (i - L) % N)
                                    for i in range(N)}
T = 2.27

filename = "configuration_profile_L_%d_T_%s.txt" % (L, T)
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

p  = 1.0 - math.exp(-2.0 / T)
nsteps = 10000
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
#    if(step % 100 == 0):
#        print "%d: M = %s" % (step / 100, M_records[-1])
E_average = sum(E_records) / float(len(E_records))
E2_average = sum(E ** 2 for E in E_records) / float(len(E_records))
c = (E2_average - E_average ** 2) / (N * T ** 2)
print 'Mean energy per spin: E/N = %s\nSpecific heat per spin: c = %s' % (E_average / N, c)

f = open(filename, 'w')
for a in S:
   f.write(str(a) + '\n')
f.close()
