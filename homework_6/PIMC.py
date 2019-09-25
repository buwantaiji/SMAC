import math, random, pylab

beta = 20.0
quartic = 1.0
cubic = - quartic

def anharmonic_potential(x, cubic, quartic):
    V = x ** 2 / 2.0 + cubic * x ** 3 + quartic * x ** 4
    return V
def anharmonic_potential_part(x, cubic, quartic):
    V_part = cubic * x ** 3 + quartic * x ** 4
    return V_part

# The (normalized) single particle density(i.e., the probability distribution p(x)) of a 1-dimension harmonic oscillator.
# Note p(x) is proportional to U(x, x, beta). Concretely, p(x) = U(x, x, beta) / Z(beta), where Z = \int dx U(x, x, beta) is the partition function.
# This function is purely used to test and check the code for the special case of harmonic potential.
def harmonic_analytic_singleparticle_density(x, beta):
    p = math.exp( - x ** 2 * math.tanh(beta / 2.0) ) \
        * math.sqrt( math.tanh(beta / 2.0) / math.pi )
    return p

# ALGORITHM 1 ###################################################################################################
# This program returns the exact free particle density matrix U_0(x, xprime, beta).
def free_density_matrix(x, xprime, beta):
    rho = math.exp( - (x - xprime) ** 2 / (2.0 * beta) ) / math.sqrt(2.0 * math.pi * beta)
    return rho

# The main program.
# For the sampling of path, we adopt a "naive" Markov Chain Monte Carlo approach.
# The initial path is arbitrarily selected.
# We propose a new path which is different from the previous path only to a single point in the path, and accept it
#   with the usual Metropolis acceptance rate.
def Metropolis_naive(M, iteration_num, x_init, iterations_per_sample, delta):
    print "RUNNING FUNCTION: Metropolis_naive"
    tau = beta / M
    x = x_init
    samples = []
    acceptance_rate = 0
    for iteration in range(iteration_num):
        k = random.randint(0, M - 1)
        k_next, k_prev = (k + 1) % M, (k - 1) % M
        x_new = x[k] + random.uniform(- delta, delta)
        weight_old =      free_density_matrix(x[k_next], x[k], tau) \
                        * free_density_matrix(x[k], x[k_prev], tau) \
                        * math.exp( - tau * anharmonic_potential(x[k], cubic, quartic) )
        weight_new =      free_density_matrix(x[k_next], x_new, tau) \
                        * free_density_matrix(x_new, x[k_prev], tau) \
                        * math.exp( - tau * anharmonic_potential(x_new, cubic, quartic) )
        if(random.uniform(0.0, 1.0) < weight_new / weight_old):
            x[k] = x_new
            acceptance_rate += 1
        if(iteration % iterations_per_sample == 0):
            samples += x
    acceptance_rate /= float(iteration_num)
    print "Acceptance_rate: %s" % acceptance_rate
    return samples
#################################################################################################################

# ALGORITHM 2 ###################################################################################################
# This function constructs a levy path with totally M slices, given start and end point x[0] = x_start and x[M] = x_end.
# The path is sampled from the probability distribution proportional to U_0(x[M], x[M-1], tau)...U_0(x[1], x[0], tau). 
#   (x[0] = x_start, x[M] = x_end is fixed; The subscript "0" stands for the free particle case.)
def free_levy_path_construction(M, tau, x_start, x_end):
    x = [x_start]
    for i in range(1, M):
        tauprime = (M - i) * tau
        mu = (tauprime * x[i - 1] + tau * x_end) \
                / (tauprime + tau)
        sigma = 1.0 / math.sqrt(1.0 / tauprime + 1.0 / tau)
        x.append( random.gauss(mu, sigma) )
    x.pop(0)
    return x

# This function constructs a total closed levy path with M slices from scratch----that is, x[0], ..., x[M-1] are all unfixed in the beginning.
# The path is sampled from the probability distribution proportional to U_0(x[0], x[M-1], tau)...U_0(x[1], x[0], tau). 
# The closeness of the path is a charateristic of the free particle partition function.
def free_levy_total_path_construction(M):
    tau = beta / M
    x_0 = random.uniform(-5.0, 5.0)
    x = [x_0] + free_levy_path_construction(M, tau, x_0, x_0)
    return x

# Algorithm2.1 ####################
# This main idea is to split the original probability distribution with anharmonic potential into two parts:
#   U(x[0], x[M-1], tau)...U(x[1], x[0], tau)  =  U_0(x[0], x[M-1], tau)...U_0(x[1], x[0], tau)
#                                               * exp( - tau * sum_{i=0}^{M-1} V(x[i]) )
# The first part correponds to the probatility distribution of free particle; while the "V(x)" in the second part stands for the whole anharmonic potential:
#   V(x) = 0.5 * x^2 + cubic * x^3 + quartic * x^4.
# In order to sample this distribution, we propose a new path according to the first part, using free levy construction;
#   Then, we accept this path with Metropolis accepatance rate corresponding to the second part.
###################################
def Metropolis_propose_total_freelevy_acceptance_fullpotential(M, iteration_num, iterations_per_sample):
    print "RUNNING FUNCTION: Metropolis_propose_TOTAL_FREElevy_acceptance_fullpotential"
    tau = beta / M
    samples = []
    x = free_levy_total_path_construction(M)
    log_weight_old =  - tau * sum(anharmonic_potential(point, cubic, quartic) for point in x)
    acceptance_rate = 0
    for iteration in range(iteration_num):
        x_new = free_levy_total_path_construction(M)
        log_weight_new =  - tau * sum(anharmonic_potential(point, cubic, quartic) for point in x_new)
        if(random.uniform(0.0, 1.0) < math.exp(log_weight_new - log_weight_old)):
            x = x_new[:]
            log_weight_old = log_weight_new
            acceptance_rate += 1
        if(iteration % iterations_per_sample == 0):
            samples += x
    acceptance_rate /= float(iteration_num)
    print "Acceptance_rate: %s" % acceptance_rate
    return samples

# Algorithm2.2 ####################
# Another version of the above algorithm, with roughly the same idea.
# Difference: In algorithm2.1, in every iteration, we propose a brand new path, with all M beads moved;
#             Now in algorithm2.2, we only propose a new path with only M_construction(1 <= M_construction <= M-1) (continuous) beads moved.
# The details of setting of propose probability and acceptance probability obeys the same idea proposed above in algorithm2.1.
# NOTE: When the parameter M_construction is suitably selected, the MC sampling can be much more effective than algorithm2.1.
#   In fact, when M_construction --> 1, this algorithm is more likely to be algorithm1. The mean acceptance rate is very HIGH, but it is likely that
#       it can't effectively can't "wander" the whole configuration space, yielding a not very efficient approach.
#   On the other hand, when M_construction --> M-1, this algorithm is more likely to be algorithm2.1.(Algorithm2.1 is actually the extreme case of M_construction = M)
#       The mean acceptance rate is very LOW, since the algorithm will then suffer from severe "overshooting" problem.
#       As a result, the MC sampling is again very inefficient, in the acceptable number of iteration time of ~ 100000 - 10000000.
#   According to the "one half" rule, the MC algorithm may achieve best performance when choosing M_construction so that the accepatance rate is roughly 50%.
###################################
def Metropolis_propose_freelevy_acceptance_fullpotential(M, M_construction, iteration_num, iterations_per_sample):
    print "RUNNING FUNCTION: Metropolis_propose_FREElevy_acceptance_fullpotential"
    tau = beta / M
    samples = []
    acceptance_rate = 0
    x = free_levy_total_path_construction(M)
    start = -1
    end = M_construction
    for iteration in range(iteration_num):
        roll = random.randint(0, M - 1)
        x = x[roll:] + x[:roll]
        log_weight_old = - tau * sum(anharmonic_potential(point, cubic, quartic) for point in x[:end])
        x_new_construction_part = free_levy_path_construction(M_construction + 1, tau, x[start], x[end])
        log_weight_new = - tau * sum(anharmonic_potential(point, cubic, quartic) for point in x_new_construction_part)
        if(random.uniform(0.0, 1.0) < math.exp(log_weight_new - log_weight_old)):
            x[:end] = x_new_construction_part[:]
            acceptance_rate += 1
        if(iteration % iterations_per_sample == 0):
            samples += x
    acceptance_rate /= float(iteration_num)
    print "Number of beads changed in an iteration: %d;  Acceptance_rate: %s" % (M_construction, acceptance_rate)
    return samples
#################################################################################################################

# ALGORITHM 3 ###################################################################################################
# This function constructs a levy path with totally M slices, given start and end point x[0] = x_start and x[M] = x_end.
# The path is sampled from the probability distribution proportional to U_{ho}(x[M], x[M-1], tau)...U_{ho}(x[1], x[0], tau). 
#   (x[0] = x_start, x[M] = x_end is fixed; The subscript "ho" stands for the harmonic oscillator case.)
def harmonic_levy_path_construction(M, tau, x_start, x_end):
    x = [x_start]
    for i in range(1, M):
        tauprime = (M - i) * tau
        a1 = 1.0 / math.tanh(tauprime) + 1.0 / math.tanh(tau)
        a2 = x[i - 1] / math.sinh(tau) + x_end / math.sinh(tauprime)
        mu = a2 / a1
        sigma = 1.0 / math.sqrt(a1)
        x.append( random.gauss(mu, sigma) )
    x.pop(0)
    return x

# This function constructs a total closed levy path with M slices from scratch----that is, x[0], ..., x[M-1] are all unfixed in the beginning.
# The path is sampled from the probability distribution proportional to U_{ho}(x[0], x[M-1], tau)...U_{ho}(x[1], x[0], tau). 
# The closeness of the path is a charateristic of the harmonic oscillator partition function.
def harmonic_levy_total_path_construction(M):
    tau = beta / M
    x_0 = random.gauss( 0.0, 1.0 / math.sqrt(2.0 * math.tanh(beta / 2.0)) )
    x = [x_0] + harmonic_levy_path_construction(M, tau, x_0, x_0)
    return x

# Algorithm3.1 ####################
# This main idea is to split the original probability distribution with anharmonic potential into two parts:
#   U(x[0], x[M-1], tau)...U(x[1], x[0], tau)  =  U_{ho}(x[0], x[M-1], tau)...U_{ho}(x[1], x[0], tau)
#                                               * exp( - tau * sum_{i=0}^{M-1} V_anharmonicpart(x[i]) )
# The first part correponds to the probatility distribution of a harmonic oscillator;
#   while the "V_anharmonicpart(x)" in the second part stands for the anharmonic part of the original potential:
#   V_anharmonicpart(x) = cubic * x^3 + quartic * x^4.
# In order to sample this distribution, we propose a new path according to the first part, using harmonic levy construction;
#   Then, we accept this path with Metropolis accepatance rate corresponding to the second part.
###################################
def Metropolis_propose_total_harmoniclevy_acceptance_anharmonicpart(M, iteration_num, iterations_per_sample):
    print "RUNNING FUNCTION: Metropolis_propose_TOTAL_HARMONIClevy_acceptance_anharmonicpart"
    tau = beta / M
    samples = []
    x = harmonic_levy_total_path_construction(M)
    weight_old = math.exp( - tau * sum(anharmonic_potential_part(point, cubic, quartic) for point in x) )
    acceptance_rate = 0
    for iteration in range(iteration_num):
        x_new = harmonic_levy_total_path_construction(M)
        weight_new = math.exp( - tau * sum(anharmonic_potential_part(point, cubic, quartic) for point in x_new) )
        if(random.uniform(0.0, 1.0) < weight_new / weight_old):
            x = x_new[:]
            weight_old = weight_new
            acceptance_rate += 1
        if(iteration % iterations_per_sample == 0):
            samples += x
    acceptance_rate /= float(iteration_num)
    print "Acceptance_rate: %s" % acceptance_rate
    return samples

# Algorithm3.2 ####################
# Another version of the above algorithm, with roughly the same idea.
# NOTE: For the motivation and insights about the differeces between these two versions of algorithm3, see the important remarks before algorithm2.2.
###################################
def Metropolis_propose_harmoniclevy_acceptance_anharmonicpart(M, M_construction, iteration_num, iterations_per_sample):
    print "RUNNING FUNCTION: Metropolis_propose_HARMONIClevy_acceptance_anharmonicpart"
    tau = beta / M
    samples = []
    x = harmonic_levy_total_path_construction(M)
    start = -1
    end = M_construction
    acceptance_rate = 0
    for iteration in range(iteration_num):
        roll = random.randint(0, M - 1)
        x = x[roll:] + x[:roll]
        weight_old = math.exp( - tau * sum(anharmonic_potential_part(point, cubic, quartic) for point in x[:end]) )
        x_new_construction_part = harmonic_levy_path_construction(M_construction + 1, tau, x[start], x[end])
        weight_new = math.exp( - tau * sum(anharmonic_potential_part(point, cubic, quartic) for point in x_new_construction_part) )
        if(random.uniform(0.0, 1.0) < weight_new / weight_old):
            x[:end] = x_new_construction_part[:]
            acceptance_rate += 1
        if(iteration % iterations_per_sample == 0):
            samples += x
    acceptance_rate /= float(iteration_num)
    print "Number of beads changed in an iteration: %d;  Acceptance_rate: %s" % (M_construction, acceptance_rate)
    return samples
#################################################################################################################

M = 100
M_construction = 20
iteration_num = 500000
x_init = [0.0] * M
iterations_per_sample = 1
delta = 1.0

#samples = Metropolis_naive(M, iteration_num, x_init, iterations_per_sample, delta)
#samples = Metropolis_propose_total_freelevy_acceptance_fullpotential(M, iteration_num, iterations_per_sample)
#samples = Metropolis_propose_freelevy_acceptance_fullpotential(M, M_construction, iteration_num, iterations_per_sample)
#samples = Metropolis_propose_total_harmoniclevy_acceptance_anharmonicpart(M, iteration_num, iterations_per_sample)
samples = Metropolis_propose_harmoniclevy_acceptance_anharmonicpart(M, M_construction, iteration_num, iterations_per_sample)

pylab.hist(samples, normed=True, bins=400, label="PIMC: algorithm 3.2 with $M_{construction}$ = %d" % M_construction)
pylab.xlim(- 3.0, 3.0)
#xgrid = [0.01 * i for i in range(-300, 301)]
#harmonic_analytic_p_x = [harmonic_analytic_singleparticle_density(x, beta) for x in xgrid]
#pylab.plot(xgrid, harmonic_analytic_p_x, label="Analytical result")
#pylab.title("(Normalized) Single particle density(i.e., probability distribution) of 1-D harmonic oscillator at inverse temperature $\\beta$\n"
#                + "$V_{ho}(x) = \\frac{1}{2} x^2$\n"
#                + "$\\beta$ = %s, $M$ = %s" % (beta, M))
pylab.title("(Normalized) Single particle density(i.e., probability distribution) of a particle in 1-D anharmonic potential at inverse temperature $\\beta$\n"
                + "$V(x) = \\frac{1}{2} x^2 + \mathrm{cubic} x^3 + \mathrm{quartic} x^4$\n"
                + "cubic = %s, quartic = %s\n" % (cubic, quartic)
                + "$\\beta$ = %s, $M$ = %s" % (beta, M))
pylab.xlabel("x")
pylab.ylabel("$p(x)$")
pylab.legend()
pylab.show()
