"""
In this experiment we propose to compare MR-systems sampled by the WE
method to control motifs sampled through the same process
"""

from weighted_ensemble_method import *
from m_r_systems import linf_norm,propose,sse,mr_system
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
    alphas = [0.5/16]*16
    mr_systems = [mr_system(alphas,use_annealing=True,scale=1,iterations=500000,
                             sse_epsilon=0.0001)[-1] for i in range(10)]
    mr_ics = map(motif_ic,mr_systems)
    mr_ginis = map(motif_gini,mr_systems)
    init_states = [random_motif_with_dirty_bits(L,num_sites)]
    ic = 12
    M = 100
    tau = 1
    timesteps = 100
    epsilon = 0.1
    bins = [-10] + myrange(0,ic,epsilon) + [ic,ic + epsilon]+ [ic + 10]
    control_motifs = concat([weighted_ensemble(mutate_motif_with_dirty_bits,
                                       motif_ic_with_dirty_bits,
                                       init_states,
                                       bins, M, tau, timesteps,
                                       final_bin_index=None,verbose=1)
                             for i in range(5)])
    control_ics = map(lambda ((m,ics),p):motif_ic(m),concat(control_motifs))
    control_ginis = map(lambda ((m,ics),p):motif_gini(m),concat(control_motifs))
    binned_controls = [[gini for ic,gini in zip(control_ics,control_ginis)
                        if abs(ic-mr_ic) < 0.1]
                       for mr_ic in mr_ics]
    
print "loaded"
