import random, math
n_trials = 400000
mean = 0.0
square_mean = 0.0
for iter in range(n_trials):
    x, y = random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0)
    Obs = 0.0
    if x**2 + y**2 < 1.0:
        Obs = 4.0
    mean += Obs
    square_mean += Obs ** 2
mean /= n_trials
square_mean /= n_trials
print "mean value: ", mean
print "standard deviation: ", math.sqrt( (square_mean - mean ** 2) * n_trials / (n_trials - 1) )
