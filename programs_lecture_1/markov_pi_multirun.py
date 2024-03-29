import random, math, pylab

def markov_pi(N, delta): 
    x, y = 1.0, 1.0
    n_hits = 0
    for i in range(N):
        del_x, del_y = random.uniform(-delta, delta), random.uniform(-delta, delta)
        if abs(x + del_x) < 1.0 and abs(y + del_y) < 1.0:
            x, y = x + del_x, y + del_y
        if x**2 + y**2 < 1.0: n_hits += 1
    return n_hits

n_runs = 1000
n_trials = 4000
delta = 0.1
ave = 0.0
std = 0.0
for run in range(n_runs):
    pi_estimate = 4.0 * markov_pi(n_trials, delta) / float(n_trials)
#    print pi_estimate
    ave += pi_estimate
    std += (pi_estimate - math.pi) ** 2
ave /= n_runs
std = math.sqrt( std / n_runs )
print "The average value = ", ave
print "The sample standard deviation = ", std
