from m_r_systems import *
from utils import *
from random import gauss

def rfreq_rseq_exp(L,G,sigma):
    matrix = [[gauss(0,sigma) for i in range(4)] for j in range(L)]
    genome = random_site(G)
    eps = [score(matrix,seq,ns=False) for seq in sliding_window(genome,L)]
    Z = sum(exp(-beta*ep) for ep in eps)
    ps = [exp(-beta*ep)/Z for ep in eps]
    r_freq = log2(G) - h(ps)
    r_seq = weighted_ic(zip(sliding_window(genome,L),ps),L)
    return r_freq,r_seq

def plot_rfreq_rseq_exp():
    L0 = 10
    G0 = 10000
    sigma0 = 1
    # Vary G:
    G_params = [(L0,G,sigma0) for G in [1000,10000,50000]]
    G_data = [[rfreq_rseq_exp(L0,G,sigma0) for i in verbose_gen(range(100))]
              for (L,G,sigma) in G_params]
    L_params = [(L,G0,sigma0) for L in [5,10,20]]
    L_data = [[rfreq_rseq_exp(L,G0,sigma0) for i in verbose_gen(range(100))]
              for (L,G,sigma) in L_params]
    sigma_params = [(L0,G0,sigma) for sigma in [0.5,1,2]]
    sigma_data = [[rfreq_rseq_exp(L0,G0,sigma) for i in verbose_gen(range(100))]
                  for (L,G,sigma) in sigma_params]
    
    plt.subplot(2,2,1)
    for i,(params,data) in enumerate(zip(G_params,G_data)):
        plt.scatter(*transpose(data),label="L=%s,G=%s,sigma=%s" % params,color="bgrycmk"[i])
    plt.legend(loc=2)
    plt.plot([0,20],[0,20])
    plt.xlabel("r_freq")
    plt.ylabel("r_seq")

    plt.subplot(2,2,2)
    for i,(params,data) in enumerate(zip(L_params,L_data)):
        plt.scatter(*transpose(data),label="L=%s,G=%s,sigma=%s" % params,color="bgr"[i])
    plt.legend(loc=2)
    plt.plot([0,20],[0,20])
    plt.xlabel("r_freq")
    plt.ylabel("r_seq")

    plt.subplot(2,2,3)
    for i,(params,data) in enumerate(zip(sigma_params,sigma_data)):
        plt.scatter(*transpose(data),label="L=%s,G=%s,sigma=%s" % params,color="bgr"[i])
    plt.legend(loc=2)
    plt.plot([0,20],[0,20])
    plt.xlabel("r_freq")
    plt.ylabel("r_seq")
    plt.show()

def weighted_ic(seqs_ps,L):
    freq_table = [[0]*4 for i in range(L)]
    for seq,p in seqs_ps:
        for i,b in enumerate(seq):
            base_index = "ACGT".index(b)
            freq_table[i][base_index] += p
    return 2*L - sum(map(h,freq_table))
    

    
print "loaded"
