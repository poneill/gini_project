"""
The purpose of this script is to explore the reasons that the greedy
sampling method produced motifs which systematically under-estimated
the true distribution of Gini coefficients over the set of all motifs
of given dimension and size.
"""

from utils import *
from collections import Counter

def beneficial_mutation_exp(trials=10000,col_length=10):
    """Count the proportion of mutations which are beneficial at each
    step in the greedy evolution of a site"""
    site = random_site(col_length)
    ic = 2 - dna_entropy(site)
    history = []
    while ic < 2:
        mutants = [mutate_site(site) for i in range(trials)]
        beneficial_mutants = filter(lambda mutant:2 - dna_entropy(mutant)>ic,
                                    mutants)
        beneficial_prop = len(beneficial_mutants)/float(trials)
        history.append((ic,beneficial_prop))
        if beneficial_mutants:
            site = beneficial_mutants[0]
            ic =  2 - dna_entropy(site)
        else:
            return history

def estimate_ben_mut_prob(site,trials=1000):
    ic = 2 - dna_entropy(site)
    mutants = [mutate_site(site) for i in range(trials)]
    beneficial_mutants = filter(lambda mutant:2 - dna_entropy(mutant)>ic,
                                mutants)
    beneficial_prop = len(beneficial_mutants)/float(trials)
    return beneficial_prop

def predict_ben_mut_prob(site):
    counts = Counter(site)
    for b in "ACGT":
        if not b in counts:
            counts[b] = 0
    n = sum(counts.values())
    l,k,j,i = sorted(counts.values())
    return (l * 1 + k*(i+j)/float(i+j+l) + j*(i/float(i+k+l)))/float(n)

def plot_beneficial_mutation_exp():
    """Plot several trajectories of beneficial mutation experiment"""
    trajectories = 10
    for _ in verbose_gen(range(trajectories)):
        history = beneficial_mutation_exp(col_length=100)
        plt.plot(*transpose(history))
    plt.show()
