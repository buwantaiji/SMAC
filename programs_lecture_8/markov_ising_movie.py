import random, math, pylab
from matplotlib.patches import FancyArrowPatch

def plot_spin(c, filename, L):
    colors = {1 : "r", -1 : "b"}
    s = 1.0 / L
    for i in range(L):
        for j in range(L):
            x, y, dy = (i + 0.5) * s, (j + 0.5) * s, 0.85 * s * c[i][j]
            arrow = FancyArrowPatch((x, y - 0.5 * dy), (x, y + 0.5 * dy),
                    fc=colors[c[i][j]], color='.2', lw=0, alpha=.8, arrowstyle="Simple" +
                    ", head_length=" + str(1.3 * 150 * s) +
                    ", head_width=" + str(1.3 * 150 * s) +
                    ", tail_width=" + str(1.3 * 40 * s))
            pylab.gca().add_patch(arrow)
    pylab.axis('scaled')
    pylab.axis([0, 1, 0, 1])
    pylab.gca().set_xticks([])
    pylab.gca().set_yticks([])
    [pylab.axhline(y=(i * s), ls='--', c='.2') for i in range(L)]
    [pylab.axvline(x=(j * s), ls='--', c='.2') for j in range(L)]
    pylab.savefig(filename)
    pylab.clf()

L = 20
N = L * L
nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
            (i // L) * L + (i - 1) % L, (i - L) % N) \
                                    for i in range(N)}
nsteps = 10000 * N
list_T = [1.0 + 0.2 * i for i in range(15)]
list_av_m = []
for T in list_T:
    S = [random.choice([1, -1]) for k in range(N)]
    M = sum(S)
    print 'T =', T
    beta = 1.0 / T
    M_tot = 0.0
    n_measures = 0
    for step in range(nsteps):
        k = random.randint(0, N - 1)
        delta_E = 2.0 * S[k] * sum(S[nn] for nn in nbr[k])
        if random.uniform(0.0, 1.0) < math.exp(-beta * delta_E):
            S[k] *= -1
            M += 2 * S[k]
#        if step % N == 0 and step > nsteps / 2:
        if True:
#            M_tot += abs(M)
            M_tot += M
            n_measures += 1
    list_av_m.append(abs(M_tot) / float(n_measures * N))
    S_lattice = [[S[i * L + j] for i in range(L)] for j in range(L)]
    plot_spin(S_lattice, "spinconf/T_%s.png" % T, L)

pylab.title('$%i\\times%i$ lattice' % (L, L))
pylab.xlabel('$T$', fontsize=16)
pylab.ylabel('$<|M|>/N$', fontsize=16)
pylab.plot(list_T, list_av_m, 'bo-', clip_on=False)
pylab.ylim(0.0, 1.0)
pylab.show()
#pylab.savefig('plot_local_av_magnetization_L%i.png' % L)
