"""
This file contains code for generating motif-recognizer (M-R) systems.

Inputs: G,n, alphas, length, sse_epsilon

Output: An M-R system whose SSE is less than epsilon.


"""
import random
from utils import product,random_site,anneal,mh,show,mutate_motif
from math import exp
from energy_matrix_analysis import score,beta

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
    
def propose(matrix,motif):
    """Return a candidate (S,R) system by mutating either the motif or
    the matrix"""
    if random.random() < 0.5:
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
    Z = Zb + sum(site_propensities)
    ps = [site_prop/Z for site_prop in site_propensities]
    return sum([(p-alpha)**2 for p,alpha in zip(ps,alphas)])

def mr_system(alphas,init_system=None,G=100000.0,n=16,L=10,
              sse_epsilon=0.0001,use_annealing=False,proposal=propose,scale=1000):
    if init_system is None:
        matrix = [[0,0,0,0] for i in range(L)]
        motif = [random_site(L) for i in range(n)]
    else:
        matrix,motif = init_system
    scaled_sse = lambda(matrix,motif):(sse(matrix,motif,alphas,G,n))*scale
    algorithm = anneal if use_annealing else mh
    if use_annealing:
        return anneal(scaled_sse,
                      lambda(matrix,motif):proposal(matrix,motif),
                      (matrix,motif),
                      iterations=500000,
                      stopping_crit = sse_epsilon*scale,return_trajectory=False)
    else:
        return mh(scaled_sse,
                  lambda(matrix,motif):proposal(matrix,motif),
                  (matrix,motif),
                  iterations=500000,
                  every=100)

def mr_system_mh(alphas,G=100000.0,n=16,L=10):
    scale = 1000 #lower means less stringent
    matrix = [[0,0,0,0] for i in range(L)]
    motif = [random_site(L) for i in range(n)]
    scaled_sse = lambda matrix,motif:(sse(matrix,motif,alphas,G,n))*scale
    return mh(lambda (matrix,motif):exp(-scaled_sse(matrix,motif)),
              lambda (matrix,motif):propose(matrix,motif),
              (matrix,motif),
              iterations=1000000,
              every=1000)

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
