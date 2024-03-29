import random, math

# The distance of two particles of position a and b, subject to periodic boundary condition.
def dist(x,y):
    d_x = abs(x[0] - y[0]) % 1.0
    d_x = min(d_x, 1.0 - d_x)
    d_y = abs(x[1] - y[1]) % 1.0
    d_y = min(d_y, 1.0 - d_y)
    return  math.sqrt(d_x**2 + d_y**2)
    
def direct_disks(N, sigma):
    n_iter = 0
    condition = False
    while condition == False:
        n_iter += 1
        L = [(random.random(), random.random())]
        for k in range(1, N):
            a = (random.random(), random.random())
            min_dist = min(dist(a, b) for b in L) 
            if min_dist < 2.0 * sigma: 
                condition = False
                break
            else:
                L.append(a)
                condition = True
    # About every n_iter attempts, there is one legal state which is sampled.
    return n_iter, L

N = 16
#eta = 0.26
eta = 0.3
sigma = math.sqrt(eta / N / math.pi)
n_runs = 100
print "Space occupation eta = %f, with corresponding sigma = %f." % (eta, sigma)
print 'Note that this program might take a while!'
for run in range(n_runs):
    iterations, config =  direct_disks(N, sigma)
#    print 'run',run
#    print iterations - 1, 'tabula rasa wipe-outs before producing the following configuration'
#    print config
#    print
    print "sample %d: n_iter = %d" % (run, iterations)

