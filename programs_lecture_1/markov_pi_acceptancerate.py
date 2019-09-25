import random, math

def markov_pi(N, delta):
    x, y = 1.0, 1.0
    n_hits = 0
    n_acceptance = 0
    for i in range(N):
        del_x, del_y = random.uniform(-delta, delta), random.uniform(-delta, delta)
        if abs(x + del_x) < 1.0 and abs(y + del_y) < 1.0:
            x, y = x + del_x, y + del_y
            n_acceptance += 1
        if x**2 + y**2 < 1.0: n_hits += 1
    acceptance_rate = n_acceptance / float(N)
    return n_hits, acceptance_rate

print "delta | acceptance rate"
for delta in [0.062, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0]:
    n_hits, acceptance_rate = markov_pi(2 ** 12, delta)
    print "%.3f" % delta + " | " + "%f" % acceptance_rate
