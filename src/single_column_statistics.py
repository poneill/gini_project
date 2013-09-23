from utils import *

def random_col(n):
    return random_site(n)

def plot_dist_of_single_col_ics(n=10,trials=10000,correct=True):
    ics = [2 - entropy(random_site(n),correct=correct,alphabet_size=4)
           for trial in xrange(trials)]
    plt.hist(ics,bins=100)

print "loaded"
