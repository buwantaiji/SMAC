import random

neighbor =  [[1, 3, 0, 0], [2, 4, 0, 1], [2, 5, 1, 2],
             [4, 6, 3, 0], [5, 7, 3, 1], [5, 8, 4, 2],
             [7, 6, 6, 3], [8, 7, 6, 4], [8, 8, 7, 5]]
t_max = 500000
site = 8
t = 0
print site
sites_statistics = [0, 0, 0, 0, 0, 0, 0, 0, 0]
while t < t_max:
    t += 1
    site = neighbor[site][random.randint(0, 3)]
#    print site
    sites_statistics[site] += 1
sites_statistics = [element / float(t_max) for element in sites_statistics]
print "The statistics of nine sites:"
print sites_statistics
