from matplotlib import pyplot as plt
from utils import *
from math import log
from alpha_sweep import *
import sys
sys.path.append("generate_coevolved_motifs")
sys.path.append("generate_coevolved_motifs/estremo_files")
from motifs import *
from estremo_utils import *
from Organism import Organism

epsilon = 10**-10

def alphas(g,n):
    return [10**-i + epsilon for i in range(10) if 10**-i > n/float(g)]

rs = [(1e-06, 0.55), (1e-06, 0.3), (1e-05, 0.57), (1e-05, 0.03),
      (0.0001, 0.06), (0.0001, 0.57), (0.001, 7.94), (0.001, 7.48),
      (0.01, 10.33), (0.01, 6.68), (0.5, 13.69), (0.5, 11.68),
      (0.99, 12.32), (0.99, 9.59)] #target occupancy, corrected Rseq
                                   #for some no-overshooting runs

sums = [(1e-06, (109000.0, 296.0)), (1e-06, (109000.0, 296.0)),
        (1e-05, (109000.0, 2210000.0)), (1e-05, (109000.0, 296.0)),
        (0.0001, (109000.0, 296.0)), (0.0001, (109000.0, 296.0)),
        (0.001, (127000.0, 1920.0)), (0.001, (131000.0, 3120.0)),
        (0.01, (400000.0, 338000.0)), (0.01, (9.02e+16, 1590000000000000.0)),
        (0.5, (109000.0, 1210000000000000.0)), (0.5, (109000.0, 1020000000000000.0)),
        (0.99, (109000.0, 1270000000000000.0)), #(target occupancy, (bg_sum,fg_sum)) 
        (0.99, (731000.0, 85600000000.0))]      #for same runs


def rfreq_star1(g,n,a):
    return log2(g) + (a*log2(a/n) + (1-a)*log2((1-a)/(g-n)))

def rfreq_star2(g,n,a):
    return log2(g*a/float(n))

def rfreq_rdagger(g,n,a):
    """What is the information of an R^dagger TF system?  I.e., what
    is the difference between the prior and posterior entropies?"""
    positives = n/a
    negatives = g - n/a
    #mass_per_positive = a/positives # = a/(n/a) = a**2/n
    mass_per_positive = 1/float(positives)
    #mass_per_negative = (1 - a)/negatives
    mass_per_negative = 1/float(negatives)
    print "positives:",positives
    print "negatives:",negatives
    print "mass_per_positive:",mass_per_positive
    print "mass_per_negative:",mass_per_negative
    print "total mass:",(positives * mass_per_positive +
                         negatives * mass_per_negative)
    assert abs(positives * mass_per_positive +
               negatives * mass_per_negative - 1) < 0.01
    h_prior = log2(g)
    h_posterior = -(positives*mass_per_positive*log2(mass_per_positive) +
                    negatives*mass_per_negative*log2(mass_per_negative))
    print "h_prior:",h_prior
    print "h_posterior:",h_posterior
    return h_prior - h_posterior


def ns_range_predicted(g,n,a):
    return log(g/float(n)) + log(a/(1-a))
    
def ns_range(fg_sum,bg_sum,n,G):
    """Given the sums over the foreground and background, compute the
    range between specific and non-specific binding energies"""
    e_fg = log(n/fg_sum)/beta
    e_bg = log((G-n)/bg_sum)/beta
    return e_bg - e_fg

def a_actual(g,n,a,rfreq_f):
    """What will be the actual probablity mass assigned to a motif by
    a recognizer using rfreq_f (i.e. a function g,n,a -> bits)?  See
    Spring 2013 notebook, date 10 March 13"""
    r_star = rfreq_f(g,n,a)
    #n * (a/n)**a * ((1-a)/float(g-n))**(1-a) # for rfreq_star1
    positives = g/float(2**r_star)
    return n/float(positives)

def rfreq_star1_is_too_low_exp(g,n,a):
    # Why do we observe ESTReMo systems using more information than
    # predicted by rfreq_star1?  I claim that this is because
    # rfreq_star1 yields information values which are simply too low
    # to accomplish the given regulatory task.

    # First let's compute R*:
    r_star = rfreq_star1(g,n,a)
    print "r_star:",r_star
    # The system allegedly requires no more than r_star bits to bind
    # to n out of g sites with probability a.  But in a genome of
    # length g, a with a recognition strength (terminology?) of b bits
    # should accept g/(2**b) sites, on average, as positives.  Hence
    # the number of positives is expected to be:
    positives = g/float(2**r_star)
    print "positives:",positives
    # Therefore the probability mass given to the true motif should be:
    a_actual = n/float(positives)
    print "a_actual:",a_actual
    efficiency = a_actual/a
    print "efficiency:",efficiency
    return efficiency
    
def plot_a_actual_exp():
    g = float(100000)
    n = 16
    alphas = [10**-i - epsilon for i in range(0,8)]
    rfreq_f = rfreq_star2
    print rfreq_f
    plt.plot(*pl(lambda a:a_actual(g,n,a,rfreq_f)/a,alphas),label="R*2")
    plt.loglog()
    plt.plot([n/g,n/g],[0.001,1])
    plt.plot([1,1],[0.001,1])
    plt.plot([n/g,1],[1,1])
    plt.xlabel("Target probability (a)")
    plt.ylabel("Efficiency (a_actual/a)")
    plt.legend(loc=0)
    plt.title("An R*2 recognizer is always efficient")
    plt.show()

def rstar_freq_1_vs_2_on_data_exp():
    n = 16
    g = float(100000)
    cutoff = n/g # point at which any recognizer can achieve target occ
    #alphas = [10**-i - epsilon for i in range(0,5)]
    #relevant_rs = [(occ,r) for occ,r in rs if occ >= cutoff]
    #plt.scatter(*transpose(relevant_rs),label="Observed")
    plt.plot(*pl(lambda a:rfreq_star1(g,n,a),myrange(0.0001,1,0.01)),label=u"$R^*$")
    plt.plot(*pl(lambda a:rfreq_star2(g,n,a),myrange(0.0001,1,0.01)),label=u"$R^\u2020$")
    plt.scatter(*transpose(alpha_sweep),label="Observed")
    plt.scatter(*transpose(alpha_sweep_0))
    plt.scatter(*transpose(alpha_sweep_3))
    #plt.plot(*pl(lambda a:rfreq_star1(g,n,a),alphas),label="R_star1")
    #plt.plot(*pl(lambda a:rfreq_star2(g,n,a),alphas),label="R_star2")
    plt.xlabel("Target Occupancy")
    plt.ylabel("$R_{seq}$ (corrected)")
    #plt.loglog()
    plt.legend(loc=4)
    plt.title("$R_{seq}$ vs. Target Occupancy")
    plt.show()
print "loaded"

def gini_motif_exp():
    from alpha_sweep_motifs import texts
    orgs = map(lambda text:Organism(text,random_site(100000),overshooting=False),
               texts)
    motifs = [org.motif for org in orgs]
    replicates = 100
    controls = [[motif_gini(control_motif(motif,epsilon=0.1,verbose=False))
                 for __ in verbose_gen(xrange(replicates))]
                for motif in verbose_gen(motifs)]
    return motifs,controls

def gini_motif_real_data_exp():
    from motifs import *
    real_motifs = [getattr(Escherichia_coli,tf_name)
                   for tf_name in Escherichia_coli.tfs]
    replicates = 100
    controls = [[motif_gini(control_motif(motif,epsilon=0.1,verbose=False))
                 for __ in verbose_gen(xrange(replicates))]
                for motif in verbose_gen(real_motifs)]
    return real_motifs,controls

def gini_motif_drosophila_exp():
    from motifs import *
    real_motifs = [getattr(Drosophila_melanogaster,tf_name)
                   for tf_name in Drosophila_melanogaster.tfs]
    replicates = 100
    controls = [[motif_gini(control_motif(motif,epsilon=0.1,verbose=False))
                 for __ in verbose_gen(xrange(replicates))]
                for motif in verbose_gen(real_motifs)]
    return real_motifs,controls

def plot_gini_motif_exp(motifs,controls,descriptor):
    motif_ginis = map(motif_gini,motifs)
    indices = sorted_indices(motif_ginis)
    sorted_motif_ginis = rslice(motif_ginis,indices)
    sorted_controls = rslice(controls,indices)
    plt.boxplot(sorted_controls)
    plt.scatter(range(1,len(motifs) + 1),sorted_motif_ginis,
             label="Observed")
    plt.ylabel("Gini Coefficient")
    plt.legend(loc=2)
    title = "Gini Coefficients for %s Motifs and Random Ensembles" % descriptor
    plt.title(title)
    
    
