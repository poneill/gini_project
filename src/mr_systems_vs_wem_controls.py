"""
In this experiment we propose to compare MR-systems sampled by the WE
method to control motifs sampled through the same process
"""

from weighted_ensemble_method import *
from m_r_systems import linf_norm,propose,sse
from utils import myrange

def sample_mr_system(alphas,G=100000.0,n=16,L=10):
    M = 100
    bins = myrange(0,0.015615,10**-5) + [1]
    tau = 1
    init_matrix = lambda:[[0,0,0,0] for i in range(L)]
    init_motif = lambda:[random_site(L) for i in range(n)]
    init_states = [(init_matrix(),init_motif()) for i in range(M)]
    q = lambda (matrix,motif):propose(matrix,motif)
    f = lambda (matrix,motif):show(sse(matrix,motif,alphas,G,n))
    timesteps = 10000
    return weighted_ensemble(q,f,init_states,bins,M,tau,timesteps,verbose=1)

def main_exp():
    L = 10
    num_sites = 16
    mr_systems = [mr_system(alphas,use_annealing=True,scale=1,iterations=500000,
                             sse_epsilon=0.0001)[-1] for i in range(30)]
    mr_ics = [motif_ic(motif) for matrix,motif in mr_systems]
    wem_ensembles = [[sample_motifs_with_dirty_bits(length=10,num_sites=16,ic=ic,
                                                   epsilon=0.1)
                      for i in verbose_gen(range(30))]
                     for ic in verbose_gen(mr_ics)]
print "loaded"
