import random, math, pylab

beta = 0.1
N = 100

# The direct sampling algorithm for a system of N distinguishable, non-interacting particles in 1-dimensional harmonic trap.
# The probability distribution p(x1, ..., xN) is proportional to U_{dist}(R, R, beta) = \prod_{i=1}^N U_{sp}(xi, xi, beta).
#   (the subscript "sp" indicates that it denotes to the density matrix of a single quantum particle in the harmonic potential.)
def harmonic_distinguishable_directsampling(iteration_num):
    samples = []
    for iterations in range(iteration_num):
        x = []
        sigma = 1.0 / math.sqrt( 2.0 * math.tanh(beta / 2.0) )
        for i in range(N):
            x.append( random.gauss(0.0, sigma) )
        samples += x
    return samples

# This function returns the (normalized) single particle density p(x) = n(x) / N for a system of N distinguishable, non-interacting
#   particles in 1-dimensional harmonic trap at inverse temperature beta.
# In this case, the single particle density p(x) is actually very trivial, and is equal to the probability density
#   in the case of a single particle in the harmonic potential.
def harmonic_distinguishable_analytic_singleparticle_density(x, beta):
    sigma = 1.0 / math.sqrt( 2.0 * math.tanh(beta / 2.0) )
    p = math.exp( - 0.5 * x ** 2 / sigma ** 2 ) / ( math.sqrt(2.0 * math.pi) * sigma )
    return p

# This function returns the density matrix U_{sp}(x, xprime, beta) of a single particle in 1-dimensional harmonic potential.
# Note: the implementation set the single particle ground state energy to zero.
#   For the original case, one should multiply the return value by exp(- beta * E_0), where E_0 = 1/2 is the ground state energy.
def harmonic_densitymatrix(x, xprime, beta):
    rho = math.exp( - (x + xprime) ** 2 * math.tanh(beta / 2.0) / 4.0
                    - (x - xprime) ** 2 / math.tanh(beta / 2.0) / 4.0 ) \
          * math.exp(beta / 2.0) \
          / math.sqrt( 2.0 * math.pi * math.sinh(beta) )
    return rho

# This function returns the single particle partition function z_{sp}(beta) of a particle in 1-dimensional harmonic potential.
# Note: this implementation set the single particle ground state energy to zero.
#   For the original case, one should multiply the return value by exp(- beta * E_0), where E_0 = 1/2 is the ground state energy.
def harmonic_z(beta):
    z = 1.0 / (1.0 - math.exp(- beta))
    return z
# This function returns the partition function Z_N(beta) of a system of N non-interacting Bosons in 1-dimensional harmonic trap.
# The analytical recursion relation is used: Z_N(beta) = 1/N * \sum{k=1}^N Z_{N-k}(beta) * z_{sp}(k*beta).
# The returned list contains the partition function of corresponding systems up to N particles: [Z_0(beta), ..., Z_N(beta)].
#   (Note Z_0(beta) = 1)
def harmonic_Boson_Z(N, beta):
    Z = [1.0]
    for n in range(1, N + 1):
        Zn = sum( harmonic_z(k * beta) * Z[n - k] for k in range(1, n + 1) ) / n
        Z.append(Zn)
    return Z
# This function returns the (normalized) single particle density p(x) = n(x) / N for a system of N non-interaction Bosons
#   in 1-dimensional harmonic trap at inverse temperature beta.
# The analytical expression is used: n(x) = 1/Z_N(beta) \sum_{k=1}^N Z_{N-k}(beta) * U_{sp}(x, x, k*beta).
#   where U_{sp}(x, x, beta) is the (diagonal) density matrix of a single particle in the harmonic potential.
#   The normalization property \int n(x) dx = N can be easily checked.
def harmonic_Boson_analytic_singleparticle_density(x, beta, N):
    Z = harmonic_Boson_Z(N, beta)
    p = sum( harmonic_densitymatrix(x, x, k * beta) * Z[N - k] for k in range(1, N + 1) ) \
        / (Z[N] * N)
    return p

# This function returns the probability distribution p(k) for a specific Boson to be in a cycle of length k(1 ~ N).
def harmonic_Boson_cycle_distribution(beta, N):
    Z = harmonic_Boson_Z(N, beta)
    Pk = [ harmonic_z(k * beta) * Z[N - k] / (N * Z[N]) for k in range(1, N + 1) ]
    return Pk

# This function returns the information related to groundstate occupation of a system of N non-interacting Bosons at inverse temperature beta.
# Return values: P:                         the probability distribution p(N0) that the single particle groundstate is occupied by N0(0 ~ N) particles.
#                mean_occupation_fraction:  The statistical mean of groundstate occupation number: <N0> / N.
def harmonic_Boson_groundstate_occupation(beta, N):
    Z = harmonic_Boson_Z(N, beta)
    P = [ (Z[N - k] - Z[N - k - 1]) / Z[N] for k in range(0, N) ]
    P.append(1.0 / Z[N])
    Z_N = Z.pop()
    mean_occupation_fraction = sum(Z) / Z_N / N
    return P, mean_occupation_fraction

def harmonic_levy_total_path_construction(M, tau):
    tau_total = M * tau
    sigma = 1.0 / math.sqrt( 2.0 * math.tanh(tau_total / 2.0) )
    x_0 = random.gauss(0.0, sigma)
    x = [x_0]
    for i in range(1, M):
        tauprime = (M - i) * tau
        a1 = 1.0 / math.tanh(tauprime) + 1.0 / math.tanh(tau)
        a2 = x[i - 1] / math.sinh(tau) + x_0 / math.sinh(tauprime)
        mu = a2 / a1
        sigma = 1.0 / math.sqrt(a1)
        x.append( random.gauss(mu, sigma) )
    return x

# The main Metropolis Monte Carlo algorithm used to simulating a system of N non-interacting Bosons.
# Several sampled quantities are returned. In a single iteration(configuration update), the following quantities of the current configuration are sampled:
#   samples: the positions of all N Bosons. One can obtain from this the single particle density of the system.
#   k_cycle_samples: the length of cycle which a randomly selected particle is in.
#   pair_distance_samples: the (absolute) distance between each pairs of the N Bosons. A histogram produced from these samples is a measure of
#                          the correlation between two Bosons of the system.
def harmonic_Boson_Metropolis(iteration_num, iterations_per_sample):
    sample_num = 0
    samples = []
    k_cycle_samples = [0] * N
    pair_distance_samples = []
    positions_permutations = {}
    sigma = 1.0 / math.sqrt( 2.0 * math.tanh(beta / 2.0) )
    for i in range(N):
        position = random.gauss(0.0, sigma)
        positions_permutations[position] = position

    for iteration in range(iteration_num):
        cycle = []
        cycle_start = random.choice( positions_permutations.keys() )
        cycle_iterator = cycle_start
        while(True):
            cycle.append(cycle_iterator)
            cycle_iterator = positions_permutations.pop(cycle_iterator)
            if(cycle_iterator == cycle_start):
                break
        cycle_length = len(cycle)
        new_cycle = harmonic_levy_total_path_construction(cycle_length, beta)
        new_cycle.append(new_cycle[0])
        for i in range(cycle_length):
            positions_permutations[ new_cycle[i] ] = new_cycle[i + 1]

        position1 = random.choice( positions_permutations.keys() )
        position1_permutation = positions_permutations.pop(position1)
        position2 = random.choice( positions_permutations.keys() )
        position2_permutation = positions_permutations.pop(position2)
        weight_old =   harmonic_densitymatrix(position1_permutation, position1, beta) \
                     * harmonic_densitymatrix(position2_permutation, position2, beta)
        weight_new =   harmonic_densitymatrix(position2_permutation, position1, beta) \
                     * harmonic_densitymatrix(position1_permutation, position2, beta)
        if(random.uniform(0.0, 1.0) < weight_new / weight_old):
            positions_permutations[position1] = position2_permutation
            positions_permutations[position2] = position1_permutation
        else:
            positions_permutations[position1] = position1_permutation
            positions_permutations[position2] = position2_permutation
        if(iteration % iterations_per_sample == 0):
            sample_num += 1
            positions = positions_permutations.keys()
            samples += positions
            k_cycle_samples[cycle_length - 1] += 1
            if(iteration % N == 0):
                pair_distance_samples += [abs(positions[i] - positions[j]) for i in range(N) for j in range(i + 1, N)]
    k_cycle_samples = [k_cycle_num / float(sample_num) for k_cycle_num in k_cycle_samples]
    return samples, k_cycle_samples, pair_distance_samples

iteration_num = 100000
iterations_per_sample = 1
#samples = harmonic_distinguishable_directsampling(iteration_num)
samples, k_cycle_samples, pair_distance_samples= harmonic_Boson_Metropolis(iteration_num, iterations_per_sample)

pylab.hist(samples, normed=True, bins=400, label="Metropolis sampling")
xgrid = [0.01 * i for i in range(-500, 501)]
distinguishable_analytic_p = [harmonic_distinguishable_analytic_singleparticle_density(x, beta) for x in xgrid]
pylab.plot(xgrid, distinguishable_analytic_p, label="Analytic result: distinguishable particles")
Boson_analytic_p = [harmonic_Boson_analytic_singleparticle_density(x, beta, N) for x in xgrid]
pylab.plot(xgrid, Boson_analytic_p, label="Analytic result: Bosons")
pylab.xlabel("x")
pylab.ylabel("$p(x) = \\frac{n(x)}{N}$")
pylab.title("Normalized single particle density of a system of $N$ non-interacting Bosons in 1-dimensional harmonic trap\n"
            + "$N$ = %d, $\\beta$ = %s" % (N, beta))
pylab.legend()
pylab.show()

################################################################
# The code below computes the probability distribution for the distance of two arbitrary particles in the system of N non-interacting Bosons
#   in 1-dimensional trap.
# Also on the same figure, the analytical result for the case of N non-interacting distinguishable particles is shown.
# One can discover the curious "Boson bunching" effect by comparing the two results.
pylab.hist(pair_distance_samples, normed=True, bins=1000, label="Metropolis sampling: Bosons")
sigma = 1.0 / math.sqrt( 2.0 * math.tanh(beta / 2.0) )
sigmaprime = math.sqrt(2) * sigma
xgrid=[0.1 * i for i in range(181)]
pair_distance_distinguishable_analytic = [ math.exp(- x ** 2 / 2.0 / sigmaprime ** 2) * 2.0 / (math.sqrt(2.0 * math.pi) * sigmaprime) for x in xgrid ]
pylab.plot(xgrid, pair_distance_distinguishable_analytic, label="Analytic result: distinguishable particles", lw=2.5)
pylab.xlabel("$r$")
pylab.ylabel("$g(r)$")
pylab.xlim(0.0, 18.0)
pylab.title("The probability distribution for the distance of two arbitrary particles in the system of N non-interacting Bosons in 1-dimensional harmonic trap\n"
            + "$N$ = %d, $\\beta$ = %s" % (N, beta))
pylab.legend()
pylab.show()

################################################################
# The code below computes the probability distribution for a specific particle to be in a permutation cycle of length k(k = 1 ~ N)
#   at various temperature.
Tstars_sparse = [0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
for Tstar in Tstars_sparse:
    beta = 1.0 / (Tstar * N)
    pylab.plot([k / float(N) for k in range(1, N + 1)],
               [p_k * N for p_k in harmonic_Boson_cycle_distribution(beta, N)],
               lw=2.5, label="$T^*$ = %s" % Tstar)
pylab.xlabel("$\\frac{k}{N}$")
pylab.ylabel("$p(k) * N$")
pylab.title("The probability distribution for a specific particle to be in a permutation cycle of length $k$(1 ~ $N$)\n"
            + "$N$ = %d, $\\beta = \\frac{1}{T^* N}$" % N)
pylab.ylim(0.0, 2.0)
pylab.legend()
pylab.show()


# The code below computes the groundstate occupation number distribution and the mean groundstate occupation number(condensate fraction)
#   at various temperature.
Tstars = [0.01 * i for i in range(1, 101)]
fraction = []
for Tstar in Tstars:
    beta = 1.0 / (Tstar * N)
    groundstate_occupation_distribution, groundstate_occupation_mean = harmonic_Boson_groundstate_occupation(beta, N)
    fraction.append(groundstate_occupation_mean)
    if(Tstar in Tstars_sparse):
        pylab.plot([occupation_num / float(N) for occupation_num in range(N + 1)],
                   [p_occupation_num * N for p_occupation_num in groundstate_occupation_distribution],
                   lw=2.5, label="$T^*$ = %s" % Tstar)
pylab.xlabel("$\\frac{N_0}{N}$")
pylab.ylabel("$p(N_0) * N$")
pylab.title("The probability distribution that the single particle groundstate is occupied by $N_0$(0 ~ $N$) particles\n"
            + "$N$ = %d, $\\beta = \\frac{1}{T^* N}$" % N)
pylab.ylim(0.0, 10.0)
pylab.legend()
pylab.show()
print "beta = %s" % beta

pylab.plot(Tstars, fraction, "r-", lw=2.5)
pylab.xlabel("$T^*$")
pylab.ylabel("$\\frac{<N_0>}{N}$")
pylab.title("The mean groundstate occupation number(condensate fraction) in a system of N non-interacting Bosons in 1-dimensional trap at inverse temperature $\\beta$\n"
            + "$N$ = %d, $\\beta = \\frac{1}{T^* N}$" % N)
pylab.show()
