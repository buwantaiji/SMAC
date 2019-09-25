import random, math

N = 4
statistics = {}
L = range(N)
nsteps = 1000000
for step in range(nsteps):
    i = random.randint(0, N - 1)
    j = random.randint(0, N - 1)
    L[i], L[j] = L[j], L[i]
    if tuple(L) in statistics: 
        statistics[tuple(L)] += 1
    else:
        statistics[tuple(L)] = 1
#    print L
#    print range(N)
#    print

samples_num = 0
print "configuration | experimental number of samples | expected number of samples"
for item in statistics:
    print item, statistics[item], float(nsteps) / math.factorial(N)
    samples_num += statistics[item]
print "The total number of samples: %d" % samples_num
