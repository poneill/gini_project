"""
This file contains code for generating motif-recognizer (M-R) systems.

Inputs: G,n, alphas, length, sse_epsilon

Output: An M-R system whose SSE is less than epsilon.


"""
import random
from utils import product,random_site,anneal,mh,show,mutate_motif
from math import exp,log
from matplotlib import rc

beta = 1.61

def score(matrix,seq,ns=True):
    """Score a sequence with a motif."""
    base_dict = {'A':0,'C':1,'G':2,'T':3}
    #specific_binding = sum([row[base_dict[b]] for row,b in zip(matrix,seq)])
    specific_binding = 0
    for i in xrange(len(matrix)):        
        specific_binding += matrix[i][base_dict[seq[i]]]
    if ns:
        return log(exp(-beta*specific_binding) + exp(-beta*ns_binding_const))/-beta
    else:
        return specific_binding            


ns_binding=False #non-specific binding is false

mono_freqs = (0.25,)*4 # uniform mononucleotide frequencies

def mutate_matrix(matrix):
    """Mutate the matrix by perturbing one weight by a standard normal"""
    L = len(matrix)
    r_i = random.randrange(L)
    r_j = random.randrange(4)
    r = random.gauss(0,0.1)
    return [[matrix[i][j]+r*(i==r_i)*(j==r_j)
             for j in range(4)] for i in range(L)]
    
def propose(matrix,motif,motif_prob=0.5):
    """Return a candidate (S,R) system by mutating either the motif or
    the matrix"""
    if random.random() < motif_prob:
        return matrix,mutate_motif(motif)
    else:
        return mutate_matrix(matrix),motif

def propose_matrix(matrix,motif):
    return mutate_matrix(matrix),motif

def propose_motif(matrix,motif):
    return matrix,mutate_motif(motif)
    
def predict_mean_prop_ref(matrix):
    """estimate <exp(-beta*score(matrix,site))>.  Ref implementation"""
    # Taken from mean_field_energy_model.py
    return product([sum([exp(-beta*ep) * freq for ep,freq in zip(col,mono_freqs)])
                    for col in matrix])

def predict_mean_prop(matrix):
    """estimate <exp(-beta*score(matrix,site))>.  See ref implementation"""
    # Taken from mean_field_energy_model.py
    sum_acc = 0
    prod_acc = 1
    for col in matrix:
        sum_acc = 0
        for ep,freq in zip(col,mono_freqs):
            sum_acc += exp(-beta*ep) * freq
        prod_acc *= sum_acc
    return prod_acc

def Z_background(matrix,G):
    """Return the contribution to the partition function from
    background binding"""
    return predict_mean_prop(matrix)*G

def sse(matrix,motif,alphas,G,n):
    Zb = Z_background(matrix,G)
    site_propensities = [exp(-beta*score(matrix,site,ns=ns_binding))
                         for site in motif]
    Z= Zb + sum(site_propensities)
    ps = [site_prop/Z for site_prop in site_propensities]
    return sum([(p-alpha)**2 for p,alpha in zip(ps,alphas)])

def linf_norm(matrix,motif,alphas,G,n):
    Zb = Z_background(matrix,G)
    site_propensities = [exp(-beta*score(matrix,site,ns=ns_binding))
                         for site in motif]
    Z= Zb + sum(site_propensities)
    ps = [site_prop/Z for site_prop in site_propensities]
    return max([abs(p-alpha) for p,alpha in zip(ps,alphas)])

def mr_system(alphas,init_system=None,G=100000.0,n=16,L=10,
              sse_epsilon=0.00000001,use_annealing=True,scale=1000,
              iterations=10000,motif_prob=0.5,verbose=False):
    proposal = lambda matrix,motif:propose(matrix,motif,motif_prob=motif_prob)
    if init_system is None:
        matrix = [[0,0,0,0] for i in range(L)]
        motif = [random_site(L) for i in range(n)]
    else:
        matrix,motif = init_system
    if use_annealing:
        scaled_sse = lambda(matrix,motif):((sse(matrix,motif,alphas,G,n))*scale)
        return anneal(scaled_sse,
                      lambda(matrix,motif):proposal(matrix,motif),
                      (matrix,motif),
                      iterations=iterations,
                      stopping_crit = sse_epsilon*scale,verbose=verbose)
    else:
        scaled_sse = lambda(matrix,motif):exp((sse(matrix,motif,alphas,G,n))*-scale)
        return mh(scaled_sse,
                  lambda(matrix,motif):proposal(matrix,motif),
                  (matrix,motif),
                  iterations=iterations,
                  every=100,verbose=True)

def mr_system_mh(alphas,G=100000.0,n=16,L=10):
    scale = 10000 #lower means less stringent
    matrix = [[0,0,0,0] for i in range(L)]
    motif = [random_site(L) for i in range(n)]
    scaled_sse = lambda matrix,motif:(sse(matrix,motif,alphas,G,n))*scale
    return mh(lambda (matrix,motif):exp(-scaled_sse(matrix,motif)),
              lambda (matrix,motif):propose(matrix,motif),
              (matrix,motif),
              iterations=100000,
              every=1000,verbose=True)


def mr_system_sa(alphas,init_system=None,G=100000.0,n=16,L=10,
              sse_epsilon=0.0001,proposal=propose,scale=1000,
              iterations=10000,return_trajectory=False):
    if init_system is None:
        matrix = [[0,0,0,0] for i in range(L)]
        motif = [random_site(L) for i in range(n)]
    else:
        matrix,motif = init_system
    scaled_sse = lambda(matrix,motif):sse(matrix,motif,alphas,G,n)*scale
    return anneal(scaled_sse,
                  lambda(matrix,motif):proposal(matrix,motif),
                  (matrix,motif),
                  iterations=iterations,
                  stopping_crit = sse_epsilon*scale,
                  return_trajectory=return_trajectory)


def alpha_ladder_exp():
    total_alphas = [1.0/2**k for k in myrange(0,10,.5)]
    alphas = [total_alpha/n for total_alpha in total_alphas]
    relative_tolerance = 0.02 # we require abs((p-a)/a) < relative_tolerance
    # follows from previous line that SSE should be less than:
    tolerances = [n*alpha**2*relative_tolerance**2 for alpha in alphas]
    xss = [[mr_system([alpha]*n,sse_epsilon=tolerance)
            for i in verbose_gen(range(10))]
           for alpha,tolerance in zip(alphas,tolerances)]
    return xss

def hold_motif_constant_exp():
    """In this experiment, we fix motifs of varying information
    contents and explore corresponding matrices via
    Metropolis-Hastings sampling."""
    ics = range(20)
    alphas = [.5/16]*16
    def evolve_trajectory(ic):
        return mh(lambda (matrix,motif):exp(-sse(matrix,motif,alphas,G,n)),
                  lambda (matrix,motif):(mutate_matrix(matrix),motif),
                  (matrix,motif),
                  iterations=10000)
    return [evolve_trajectory(ic) for ic in verbose_gen(ics)]
    
print "loaded mr_systems"

def mh_motif(n,w,desired_ic,epsilon,scale=10,iterations=10000):
    """Find a motif satisfying desired_ic +/- epsilon by mh sampling"""
    motif = [random_site(w) for i in range(n) ]
    f = lambda m:exp(-abs(desired_ic-motif_ic(m))*scale)
    proposal = mutate_motif
    return mh(f,proposal,motif,iterations=iterations)

def mh_control_motif(motif,epsilon=0.1):
    n = len(motif)
    w = len(motif[1])
    desired_ic = motif_ic(motif)
    return mh_motif(n,w,desired_ic,epsilon,scale=20)[-1]
    
def mh_motif_phase_transition_fig(filename=None):
    rc('text',usetex=True)
    c = 5
    for s in range(1,51):
        print s
        plt.scatter(s,mean(map(lambda m:(motif_ic(m)-c)**2,mh_motif(16,10,c,.1,scale=s,iterations=5000))))
    plt.semilogy()
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$<(IC(M)-c)^2>$ (bits)$^2$")
    maybesave(filename)
    
def sa_motifs_vs_controls_gini_comparison():
    sa_motifs = [mr_system_sa([.5/16 for i in range(16)],scale=1,
                           iterations=1000000,return_trajectory=False)[1]
              for i in verbose_gen(range(100))][:20]
    motif_ginis = map(motif_gini,sa_motifs)
    replicates = 100
    control_ginis = [[motif_gini(control_motif(motif,epsilon=0.1))
                      for __ in verbose_gen(xrange(replicates))]
                     for motif in verbose_gen(sa_motifs[:20])]
