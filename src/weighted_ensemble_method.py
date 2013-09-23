"""
The purpose of this script is to attempt to sample uniformly from the
set of all motifs having IC of c +/- epsilon.  We do this through the
weighted ensemble method.
"""
from random import shuffle,random,choice
from utils import *

def total_prob(binned_states):
    return sum(p for bs in binned_states for (state,p) in bs)

def weighted_ensemble(q, f, init_states, bins, M, tau, timesteps):
    Q = lambda (state,p):(q(state),p)
    F = lambda (state,p):(f(state),p)
    n = float(len(init_states))
    binned_states = [[(state,1/n) for state in init_states if b0 <= f(state) < b1]
              for b0,b1 in pairs(bins)]
    #print "total prob:",total_prob(binned_states)
    t = 0
    #history = []
    while t < timesteps:
        next_tau = t + tau
        while t < next_tau:
            # Run the dynamics
            binned_states = mmap(Q,binned_states)
            t += 1
        # resort trajectories to correct bins
        # inefficient!  quadratic!
        # binned_states = [[(state,p) for binned_state in binned_states
        #                   for (state,p) in binned_state
        #                   if b0 <= f(state) < b1]
        #                  for b0,b1 in pairs(bins)]
        new_binned_states = [[] for _ in pairs(bins)]
        for bs in binned_states:
            for state,p in bs:
                f_state = f(state)
                for i,(b0,b1) in enumerate(pairs(bins)):
                    if b0 <= f_state <b1:
                        index = i
                        break
                new_binned_states[i].append((state,p))
        binned_states = new_binned_states
        # trajectories are now resorted into correct bins
        for bs_index in range(len(binned_states)):
            binned_state = binned_states[bs_index]
            # If bin must be upsampled...
            if 0 < len(binned_state) < M: 
                deficiency = M - len(binned_state)
                p_before = sum([p for state,p in binned_states[bs_index]])
                available_states = binned_states[bs_index]
                total_p = sum([p for m,p in available_states])
                ps = [p/total_p for m,p in available_states]
                new_states = [inverse_cdf_sample(available_states,ps)
                              for i in range(deficiency)]
                binned_states[bs_index].extend(new_states)
                p_after = sum([p for state,p in binned_states[bs_index]])
                binned_states[bs_index] = [(state,p*p_before/p_after)
                                           for (state,p) in binned_states[bs_index]]
                assert abs(p_before - p_after) < 10**-50
            #bin must be downsampled...
            elif len(binned_state) > M:
                surplus = len(binned_state) - M
                p_before = sum([p for state,p in binned_states[bs_index]])
                for i in range(surplus):
                    binned_states[bs_index].pop()
                p_after = sum([p for state,p in binned_states[bs_index]])
                binned_states[bs_index] = [(state,p*p_before/p_after)
                                           for state,p in binned_states[bs_index]]
        probs = map(lambda bs:sum(p for s,p in bs),binned_states)
        print t,probs
    return binned_states

def sample_motifs(length,num_sites,ic,epsilon):
    bins = [-10] + myrange(0,ic,epsilon) + [ic,ic + epsilon]+ [ic + 10]
    tau = 1
    timesteps = 50
    results = weighted_ensemble(q, f, init_states, bins, M, tau, timesteps)
    motifs,ps = transpose(concat(results[-3:-1]))
    return inverse_cdf_sample(motifs,normalize(ps))

def wem_control_motif(motif,epsilon=0.1):
    num_sites = len(motif)
    length = len(motif[0])
    ic = motif_ic(motif)
    return sample_motifs(length,num_sites,ic,epsilon)

def test_we():
    """Check the validity of the implementation by sampling the random
    walk on the 1d integer lattice"""
    min_bin = 0
    max_bin = 100
    def q(x):
        if x == min_bin:
            return min_bin + 1
        elif x == max_bin:
            return max_bin - 1
        else:
            return x- 1 if random.random() < .5 else x + 1
    f = lambda x:x
    init_states = [0]
    bins = range(0,max_bin,max_bin/10) + [100]
    M = 1000
    tau = 1000
    timesteps = tau * 100
    results = weighted_ensemble(q, f, init_states, bins, M, tau, timesteps)
    return results

def ic_we():
    bins = range(-1,3) + myrange(3,12,0.1) + [50]
    hist = weighted_ensemble(mutate_motif,
                             motif_ic,
                             [random_motif(10,16) for i in range(100)],
                             bins=bins,
                             M=1000,
                             tau=1,
                             timesteps=500)
    return hist

print "loaded"
