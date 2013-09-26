from utils import *
from weighted_ensemble_method import *
import sys
sys.path.append("../results/coevolved_motifs/estremo_motifs")
from estremo_motifs import estremo_motifs
import cPickle

def estremo_motif_exp(filename=None):
    # Since the estremo motifs are all of equivalent dimension, we can
    # use the entire distribution of binned states
    length = 10
    num_sites = 16
    estremo_ics = map(motif_ic,estremo_motifs)
    ic = max(estremo_ics)
    M = 100
    timesteps = 100
    tau = 1
    epsilon = 0.1
    init_states = [random_motif_with_dirty_bits(length,num_sites)
                   for i in range(M)]
    bins = [-10] + myrange(0,ic,epsilon) + [ic,ic + epsilon]+ [ic + 10]
    controls = transpose(concat([concat(weighted_ensemble(mutate_motif_with_dirty_bits,
                                                          motif_ic_with_dirty_bits,
                                                          init_states,
                                                          bins, M, tau, timesteps))
                                 for i in range(10)]))[0]
    control_ics = map(motif_ic,controls)
    matched_controls = [[control_motif for control_motif,control_ic
                         in zip(controls,control_ics) if abs(control_ic - ic) < epsilon]
                        for ic in estremo_ics]
    matched_ginis = mmap(motif_gini,matched_controls)
    estremo_ginis = map(motif_gini,estremo_motifs)
    indices = sorted_indices(estremo_ics)
    plt.boxplot(rslice(matched_ginis,indices))
    plt.scatter(range(1,len(estremo_motifs) + 1),rslice(estremo_ginis,indices),label="observed")
    plt.ylabel("Gini Coefficient")
    plt.legend(loc=0)
    plt.title("Gini Coefficients for ESTReMo Motifs vs. Random Ensembles")
    maybesave(filename)
    with open("../results/estremo_controls.pkl",'w') as f:
        cPickle.dump(controls,f)

def ecoli_motif_exp():
    ecoli_motifs = [getattr(Escherichia_coli,tf) for tf in Escherichia_coli.tfs]
    dim = lambda m:(len(m),len(m[0]))
    dims = map(dim,ecoli_motifs)
    cutoff = 10
    sorted_js = sorted_indices(dims)
    ecoli_controls = [match_motif(motif,inv_cdf=False)
                      for motif in rslice(ecoli_motifs,sorted_js)[:cutoff]]
    ecoli_ginis = map(motif_gini,ecoli_motifs)
    ecoli_control_ginis = mmap(motif_gini,ecoli_controls)
    plt.boxplot(ecoli_control_ginis)
    plt.scatter(range(1,cutoff+1),rslice(ecoli_ginis,sorted_js)[:10])
    with open("../results/ecoli_controls.pkl",'w') as f:
        cPickle.dump(econli_controls,f)
