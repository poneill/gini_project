"""
We wish to explore the relationship between copy number and various
energy matrix statistics, such as disorder (i.e. standard deviation of
score).
"""

from utils import *
from parse_energy_matrices import parse_energy_matrices
from scipy.stats import wilcoxon
#energy_matrices = parse_energy_matrices("energy_matrices.txt")

beta = 1.61
n = 16
G = 100000.0
alpha = 0.5
L = 10
ns_binding_const = -8
base_dict = {b:i for i,b in enumerate("ACGT")}

def score_ref(matrix,seq,ns=True):
    """Score a sequence with a motif.  This is a reference
    implementation.  See score for production implementation"""
    specific_binding = sum([row[base_dict[b]] for row,b in zip(matrix,seq)])
    if ns:
        return log(exp(-beta*specific_binding) + exp(-beta*ns_binding_const))/-beta
    else:
        return specific_binding

def score(matrix,seq,ns=True):
    """Score a sequence with a motif."""
    #specific_binding = sum([row[base_dict[b]] for row,b in zip(matrix,seq)])
    specific_binding = 0
    for i in xrange(len(matrix)):        
        specific_binding += matrix[i][base_dict[seq[i]]]
    if ns:
        return log(exp(-beta*specific_binding) + exp(-beta*ns_binding_const))/-beta
    else:
        return specific_binding            
    
def matrix_mean(matrix):
    """Return the mean score for the energy matrix"""
    return sum(map(mean,matrix))

def matrix_variance(matrix):
    """Return the variance of the scores for the energy matrix"""
    return sum(map(lambda row:variance(row,correct=False),matrix))

def matrix_sd(matrix):
    return sqrt(matrix_variance(matrix))

def specific_binding_fraction(matrix,n=10000):
    """What fraction of the time does the tf bind specifically (i.e. < -8kbt)
    to a random site?"""
    return mean([score(matrix,random_site(10)) < -8 for i in xrange(n)])

def predict_mean_prop(matrix,ns=True):
    """estimate <exp(-beta*score(matrix,site))>
    ns: non-specific binding"""
    return (product([mean([exp(-beta*ep) for ep in col]) for col in matrix]) +
            (exp(-beta*(-8)) if ns else 0)) # kbT
            

def predict_variance_prop(matrix):
    """estimate Var(exp(-beta*score(matrix,site)))"""
    # See: The Variance of the Product of K Random Variables
    # Leo A. Goodman
    # Page 55 of 54-60

    # However, first line of equation 1 is /wrong/.  Should be:
    # product(e_of_sqs) - product(esqs), not sum...

    # expectation of square
    e_of_sqs = [mean([exp(-beta*ep)**2 for ep in col]) for col in matrix]
    # square of expectation
    esqs = [mean([exp(-beta*ep) for ep in col])**2 for col in matrix]
    return product(e_of_sqs) - product(esqs)

def predict_z(matrix,num_sites,ns=True):
    return predict_mean_prop(matrix,ns=ns) * num_sites

def predict_z_variance(matrix,num_sites):
    return predict_variance_prop(matrix) * num_sites

def mean_variance_plot(filename=None):
    means = map(matrix_mean,energy_matrices)
    variances = map(matrix_variance,energy_matrices)
    plt.scatter(ks,means,label="Mean")
    plt.scatter(ks,variances,label="Variance",color='g')
    mean_regression = lambda x:poly1d(polyfit(map(log,ks),means,1))(log(x))
    variance_regression = lambda x:poly1d(polyfit(map(log,ks),variances,1))(log(x))
    plt.plot(*pl(mean_regression,map(iota,range(1,65))))
    plt.plot(*pl(variance_regression,map(iota,range(1,65))))
    plt.semilogx()
    plt.xlabel("Copy Number")
    plt.ylabel("kBT,(kBT)^2")
    plt.legend(loc=0)
    if filename:
        plt.savefig(filename,dpi=400)
    else:
        plt.show()

def specific_binding_fraction_plot(filename=None):
    binding_fractions = map(specific_binding_fraction,verbose_gen(energy_matrices))
    plt.scatter(ks,binding_fractions)
    plt.xlabel("Copy Number")
    plt.ylabel("Specific binding fraction")
    plt.loglog()
    if filename:
        plt.savefig(filename,dpi=400)
    else:
        plt.show()
    
def max_16_ic(matrix,n=100000):
    """
    Compute motif_ic of top 16 of n random sites
    """
    width = len(matrix)
    sites = [random_site(width) for i in xrange(n)]
    scores = map(lambda site:score(matrix,site),sites)
    top16 = map(first,sorted(zip(sites,scores),key=lambda(site,score):score)[:16])
    return motif_ic(top16)

def matrix_as_psfm(matrix):
    """
    convert energy to psfm, assuming uniform probabilities
    DEPRECATED: INCOMPLETE
    """
    return [[2**(ep-2) for ep in row] for row in matrix]

def predict_site_energy(matrix,n,G,alpha):
    """See blue notebook: 6/27/13"""
    #constant which does not depend on sites, matrix
    C = 1/beta * log(n/G*(1-alpha)/alpha)
    Zb = predict_z(matrix,G)
    omega = Zb/G
    return C - (1/beta * log(omega))

def predict_zf(matrix,n,G,alpha):
    """Predict sum_{i=1}^n exp(-\beta*E(s))"""
    ep_f = predict_site_energy(matrix,n,G,alpha)
    return n*exp(-beta*ep_f)

C = 1/beta * (log(n/G*(1-alpha)/alpha) + L * log(4))
def site_error(matrix,site):
    """Compute error for site, given matrix"""
    return (score(matrix,site,ns=False)
            - C
            + 1/beta * sum([log(sum([exp(-beta*matrix[i][base_dict[b]])
                                   for b in "ACGT"]))
                           for i in range(L)]))

def site_error_optimized(matrix,site):
    """Compute error for site, given matrix"""
    return score(matrix,site,ns=False) - C

def sse(matrix,motif):
    """Compute sum of squared error for matrix and motif"""
    return sum([site_error(matrix,site)**2
                for site in motif])

def sse_optimized(matrix,motif):
    """Compute sum of squared error for matrix and motif"""
    #Hoisted computation of K out of site_error
    K = 1/beta * sum([log(sum([exp(-beta*matrix[i][base_dict[b]])
                                   for b in "ACGT"]))
                           for i in range(L)])
    return sum([(site_error_optimized(matrix,site)+K)**2
                for site in motif])

def sse_experiment(motif):
    """Given a collection of sites, can we find a corresponding energy
    matrix by gradient descent on sum of squared errors?"""
    L = len(matrix)
    n = len(motif)
    G = 100000.0
    alpha = 0.9
    
    def partial_site_error(matrix,site,i,b):
        return 2*site_error(matrix,site)*(int(site[i] == b)
                                          - (exp(-beta*matrix[i][base_dict[b]])
                                             /sum([exp(-beta*matrix[i][base_dict[c]])
                                                   for c in "ACGT"])))
    
    
    def partial_sse(matrix,i,b):
        return sum([2*site_error(matrix,site)*partial_site_error(matrix,site,i,b)
                    for site in motif])

    def jacobian(matrix):
        return [[partial_sse(matrix,i,b) for b in "ACGT"] for i in range(L)]

    def grad_desc(matrix,max_its=1000):
        step = 0.0001
        tolerance = 0.001
        current_sse = sse(matrix)
        print current_sse
        its = 0 # iterations
        sses = []
        while current_sse > tolerance and its < max_its:
            j = jacobian(matrix)
            #print j
            matrix = matrix_add(matrix,matrix_scalar_mult(-step,j))
            current_sse = sse(matrix)
            print its,current_sse #,[score(matrix,site,ns=False) for site in motif]
            #print matrix
            its += 1
            sses.append(current_sse)
        return matrix,sses

    print matrix,motif
    return grad_desc([[0]*4 for i in range(L)])

def mutate_matrix(matrix):
    """Mutate the matrix by perturbing one weight by a standard normal"""
    L = len(matrix)
    r_i = random.randrange(L)
    r_j = random.randrange(4)
    r = random.gauss(0,1)
    return [[matrix[i][j]+r*(i==r_i)*(j==r_j)
             for j in range(4)] for i in range(L)]
    
def propose(matrix,motif):
    """Return a candidate (S,R) system by mutating either the motif or
    the matrix"""
    if random.random() < 0.5:
        return matrix,mutate_motif(motif)
    else:
        return mutate_matrix(matrix),motif

def mh_experiment(text="",filename=None):
    """Metropolis-Hastings sampling for SSE of (S,R) systems"""
    motif = [random_site(L) for i in range(n)]
    matrix = [[0,0,0,0] for i in range(L)]
    xs = mh(lambda(matrix,motif):exp(-sse_optimized(matrix,motif)), 
            lambda(matrix,motif):propose(matrix,motif),
            (matrix,motif),
            verbose=True,
            iterations=50000)
    sses = [sse_optimized(matrix,motif)
            for (matrix,motif) in verbose_gen(xs,modulus=1000)] 
    ics = [motif_ic(motif) for (matrix,motif) in verbose_gen(xs,modulus=1000)]
    ginis = [motif_gini(motif) for (matrix,motif) in verbose_gen(xs,modulus=1000)]
    plt.scatter(ics,sses)
    plt.xlabel("Motif Information Content (bits)")
    plt.ylabel("Sum Squared Error")
    plt.title("Motif IC vs. Sum Squared error: %s " % text)
    maybesave(filename)
    
def ic_vs_gini_scatterplot_exp(trials=100,stopping_crit=1,filename=None):
    # redo this properly!
    matrix = [[0,0,0,0] for i in range(L)]
    motif = [random_site(L) for i in range(n)]
    # xs = mh(lambda(matrix,motif):exp(-sse_optimized(matrix,motif)),
    #               lambda(matrix,motif):propose(matrix,motif),
    #               (matrix,motif),
    #               verbose=True,
    #               iterations=50000,every=100)
    scale = 0.01 #use this to prevent overflows in anneal
    scaled_sse = lambda(matrix,motif):sse_optimized(matrix,motif)*scale
    annealed_system = lambda: anneal(scaled_sse,
                                     lambda(matrix,motif):propose(matrix,motif),
                                     (matrix,motif),
                                     iterations=100000,
                                     stopping_crit = stopping_crit*scale)
    xs = [annealed_system() for i in verbose_gen(xrange(trials))]
    ginis = [motif_gini(motif)
                   for (matrix,motif) in verbose_gen(xs,modulus=1000)]
    ics = [motif_ic(motif) for (matrix,motif) in verbose_gen(xs,
                                                             modulus=1000)]
    
    sa_motifs = [sa_motif_with_desired_ic(ic,0.1,16,10)
                         for ic in verbose_gen(ics)]
    greedy_motifs = [generate_greedy_motif_with_ic(ic,0.1,16,10)
                             for ic in verbose_gen(ics)]
    
    sa_ics = [motif_ic(motif) for motif in sa_motifs]
    sa_ginis = [motif_gini(motif) for motif in sa_motifs]

    greedy_ics = [motif_ic(motif) for motif in greedy_motifs]
    greedy_ginis = [motif_gini(motif) for motif in greedy_motifs]
    
    plt.scatter(ics,ginis,color='b',label="Systems")
    plt.scatter(sa_ics,sa_ginis,color='r',label="MCMC")
    plt.scatter(greedy_ics,greedy_ginis,color='g',label="Greedy")
    print "Systems vs. SA motifs:",wilcoxon(ginis,sa_ginis)
    print "Systems vs. greedy motifs:",wilcoxon(ginis,greedy_ginis)
    print "Greedy vs. sa motifs:",wilcoxon(sa_ginis,greedy_ginis)
    # Systems vs. SA motifs: (75677.0, 2.1175931461637028e-81)
    # Systems vs. greedy motifs: (78734.0, 1.2196302463732947e-78)
    # Greedy vs. sa motifs: (247043.0, 0.72555385752004398)

    maybesave(filename)
    
print "loaded energy_matrix_analysis"

def mr_pairs_have_less_mi_exp(filename=None):
    """Motifs evolved as M-R pairs have less MI than random motifs
    with the same IC"""
    trials = 500
    matrix = [[0,0,0,0] for i in range(L)]
    motif = [random_site(L) for i in range(n)]
    scale = 0.01 #use this to prevent overflows in anneal
    scaled_sse = lambda(matrix,motif):sse_optimized(matrix,motif)*scale
    annealed_system = lambda :anneal(scaled_sse,
                                     lambda(matrix,motif):propose(matrix,motif),
                                     (matrix,motif),
                                     verbose=True,
                                     iterations=100000,
                                     stopping_crit = 0.1*scale)
    systems = [annealed_system() for i in xrange(500)]
    motifs = map(second,systems)
    ics = map(motif_ic,motifs)
    control_motifs = [sa_motif_with_desired_ic(ic,0.1,n,L) for ic in verbose_gen(ics)]
    mis = map(total_motif_mi,motifs)
    control_mis = map(total_motif_mi,control_motifs)
    plt.scatter(mis,control_mis)
    plt.xlabel("M-R System Mutual Information (bits)")
    plt.ylabel("Annealed Motif Mutual Information (bits)")
    plt.plot([0,5],[0,5])
    maybesave(filename)
    #mannwhitneyu(mis,control_mis) -> (47673.0, 1.2864021557444156e-64)
    return mis,control_mis

