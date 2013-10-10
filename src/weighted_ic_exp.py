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

def random_matrix(L,sigma):
    return [[random.gauss(0,sigma) for i in range(4)] for j in range(L)]

def matrix_variance(matrix,centered=True):
    return sum(mean([x**2 for x in row]) - centered*mean(row)**2
               for row in matrix)

def matrix_variance_check(matrix):
    return sum(3/16.0*sum(matrix[i][b]**2 for b in range(4)) -
               1/16.0*sum(matrix[i][b]*matrix[i][c] for b in range(4)
                          for c in range(4) if not b == c)
               for i in range(len(matrix)))

def variance_exp(independent_matrices=False):
    trials = 100
    G = 10000
    variance_data = []
    for trial in verbose_gen(xrange(trials)):
        L = random.randrange(1,100)
        sigma = random.random() * 10
        genome = random_site(G)
        if not independent_matrices:
            matrix = random_matrix(L,sigma)
            eps = [score(matrix,seq,ns=False) for seq in sliding_window(genome,L)]
        else:
            eps = [score(random_matrix(L,sigma),seq,ns=False)
                   for seq in sliding_window(genome,L)]
        variance_data.append((variance(eps),L,sigma))
    return variance_data
        
def plot_var_data(var_data):
    plt.scatter(*transpose([(3/4.0*L*sigma**2,obs) for obs,L,sigma in var_data]))
    max_obs = max(map(first,var_data))
    plt.plot([0,max_obs],[0,max_obs])
    plt.plot([0,max_obs],[0,4/3.0*max_obs])

def variance_underprediction_caused_by_using_one_matrix():
    dependent_data = variance_exp()
    independent_data = variance_exp(independent_matrices=True)
    plot_var_data(dependent_data)
    plot_var_data(independent_data)
    plt.show()

def expectation_of_Z_exp():
    L = 10
    G = 10000
    sigma = 2
    matrix = random_matrix(L,sigma)
    genome = random_site(G+L-1)
    eps = [score(matrix,seq,ns=False)
           for seq in sliding_window(genome,L)]
    print "normality test:",normaltest(eps)
    print "eps variance:",variance(eps)
    print "predicted eps variance:",(3/4.0*L*sigma**2)
    assert(len(eps) == G)
    print "eps variance w/ beta:",variance(map(lambda x:x*-beta,eps))
    Sigma_sq = 3/4.0*L*sigma**2*beta**2
    print "predicted eps variance w/ beta:",(Sigma_sq)
    Z = sum(exp(-beta*ep) for ep in eps)
    expected_Z = G*exp((Sigma_sq)/2)
    return Z,expected_Z

def lognormal_test():
    mu = 1
    sigma = 2
    c = 0.1
    xs = [random.gauss(mu,sigma) for _ in xrange(10000)]
    ys = [exp(c*x) for x in xs]
    y_bar = mean(ys)
    y_pred = exp(c*mu+c**2*sigma**2/2)
    return y_bar,y_pred

print "loaded"
