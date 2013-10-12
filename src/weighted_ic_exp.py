from m_r_systems import *
from utils import *
from random import gauss
from scipy.stats import norm

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

def matrix_mean(matrix):
    return sum(mean([x for x in row])
               for row in matrix)

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

def mean_eps_exp(G,L,sigma):
    """
    Conclusion: <ep> ~= min(eps) for sigma >>1
    """
    matrix = [[gauss(0,sigma) for i in range(4)] for j in range(L)]
    genome = random_site(G-L+1)
    eps = [score(matrix,seq,ns=False) for seq in sliding_window(genome,L)]
    Z = sum(exp(-beta*ep) for ep in eps)
    ps = [exp(-beta*ep)/Z for ep in eps]
    mean_energy = sum(p*ep for p,ep in zip(ps,eps))
    return min(eps),mean_energy

def mean_eps_exp2(mu,sigma,n):
    eps = [random.gauss(mu,sigma) for i in range(n)]
    Z = sum(exp(-beta*ep) for ep in eps)
    ps = [exp(-beta*ep)/Z for ep in eps]
    mean_ep = sum(ep*p for ep,p in zip(eps,ps))
    Z_pred = n*exp(-beta*mu+beta**2*sigma**2/2.0)
    db = 0.01
    dZ_pred = -(n*exp(-(beta + db)*mu+(beta+db)**2*sigma**2/2.0) - Z_pred)/db
    print "Z:",Z,"Z_pred:",Z_pred
    pred_mean_ep = (-n*exp(-beta*mu+(beta*sigma)**2/2.0)*(-mu+beta*sigma**2))/Z_pred
    print "mean ep:",mean_ep,"pred_mean_ep:",pred_mean_ep,"dZ_pred",dZ_pred
    return pred_mean_ep,mean_ep
    
    
    

def logZ(G,L,sigma):
    matrix = [[gauss(0,sigma) for i in range(4)] for j in range(L)]
    genome = random_site(G-L+1)
    eps = [score(matrix,seq,ns=False) for seq in sliding_window(genome,L)]
    matrix_mu = matrix_mean(matrix)
    matrix_sigmasq = matrix_variance(matrix)
    Z = sum(exp(-beta*ep) for ep in eps)
    return log(Z)

def logZ_exp(trials=100):
    # Fenton-Wilkinson method: 
    # ep~N(0,3/4*L*sigma**2), so assume
    # logZ ~ G*N(0,3/4*L*sigma**2)
    G = 10000
    L = 10
    sigma = 1
    matrix = [[gauss(0,sigma) for i in range(4)] for j in range(L)]
    matrix_mu = matrix_mean(matrix)
    matrix_sigmasq = matrix_variance(matrix)
    Xss = [[score(matrix,seq,ns=False)
                  for seq in sliding_window(random_site(G-L+1),L)]
           for i in verbose_gen(range(trials))]
    Yss = [[exp(-beta*ep) for ep in Xs] for Xs in Xss]
    logZs = [log(sum(Ys))
             for Ys in Yss]
    #xs = [(beta*rlognorm(matrix_mu,matrix_sigmasq)) for t in range(trials)]
    #qqplot(logZs,xs)
    pred = exp((-beta)*matrix_mu+(beta**2)*matrix_sigmasq/2)
    print pred,mean(concat(Yss))
    plt.boxplot(logZs)
    plt.plot([0,1],[pred,pred])
    return logZs,pred

def rlognorm(mu,sigmasq):
    return exp(random.gauss(mu,sqrt(sigmasq)))

def dlognorm(x,mu,sigma):
    return 1/(x*sqrt(2*pi)*sigma)*exp(-(log(x)-0)**2/(2*sigma**2))

def h_logZ_exp(G,L,sigma):
    """
    Conclusion: posterior entropy (in nats) == beta*mean_energy+log(Z)
    This result (really a stat mech identity) is a key step in
    predicting post ent.  This implies that tractable predictions for
    mean_energy, logZ are crucial.
    """
    matrix = [[gauss(0,sigma) for i in range(4)] for j in range(L)]
    genome = random_site(G-L+1)
    eps = [score(matrix,seq,ns=False) for seq in sliding_window(genome,L)]
    Z = sum(exp(-beta*ep) for ep in eps)
    ps = [exp(-beta*ep)/Z for ep in eps]
    assert(abs(sum(ps) - 1) < 10**-10)
    # return entropy in nats, log(Z)
    mean_energy = sum(p*ep for p,ep in zip(ps,eps))
    H = h(ps)*log(2) #entropy in nats
    H2 = -sum(p*(-beta*ep-log(Z)) for p,ep in zip(ps,eps))
    H3 = sum(p*(beta*ep)+p*log(Z) for p,ep in zip(ps,eps))
    H4 = beta*mean_energy + sum(p*log(Z) for p in ps)
    Hfinal = beta*mean_energy + log(Z)
    print H,H2,H3,H4,Hfinal
    #return (h(ps)*log(2),beta*mean_energy+log(Z))
    return (h(ps)*log(2),beta*mean_energy+log(Z))

def predict_post_ent(G,L,sigma,matrix_stats=False):
    matrix = [[gauss(0,sigma) for i in range(4)] for j in range(L)]
    genome = random_site(G-L+1)
    eps = [score(matrix,seq,ns=False) for seq in sliding_window(genome,L)]
    Z = sum(exp(-beta*ep) for ep in eps)
    ps = [exp(-beta*ep)/Z for ep in eps]
    H = h(ps)*log(2) #entropy in nats
    matrix_mu = 0 if not matrix_stats else matrix_mean(matrix)
    matrix_sigma = 3/4.0*L*sigma**2 if not matrix_stats else sqrt(matrix_variance(matrix))
    pred_min_ep = matrix_mu + norm.ppf(1/float(G+1))*matrix_sigma
    pred_logZ = log(G) + -beta*matrix_mu + (beta*matrix_sigma)**2/2.0
    print "---"
    print "G:",G,"L:",L,"sigma:",sigma
    print "logZ:",pred_logZ,log(Z)
    mean_energy = sum(p*ep for p,ep in zip(ps,eps))
    print "energy:","mean:",mean_energy,"min:",min(eps),"pred_min:",pred_min_ep
    pred_H = beta*pred_min_ep + pred_logZ
    print "H:",pred_H,H
    return pred_H,H

    
    
def lognormal_test(mu=1,sigma=2,c=1):
    xs = [random.gauss(mu,sigma) for _ in xrange(10000)]
    ys = [exp(c*x) for x in xs]
    y_bar = mean(ys)
    y_pred = exp(c*mu+c**2*sigma**2/2)
    return y_bar,y_pred

def powersum_test_ref(mu=1,sigma=2,n=10000):
    logz = log(sum(exp(random.gauss(mu,sigma)) for _ in xrange(n)))
    logz_pred = log(n) + mu+sigma**2/2.0
    return logz_pred,logz

def powersum_test(mu=1,sigma=2,n=10000):
    acc = 0
    for _ in xrange(n):
        acc += exp(random.gauss(mu,sigma))
    logz = log(acc)
    logz_pred = log(n) + mu+sigma**2/2.0
    return logz_pred,logz

def plot_powersum_test(G=10000,trials=100):
    mu_range = 5 #[-mu_range,mu_range]
    sigma_range = 10 #[0,sigma_range]
    psums = [powersum_test(0,#random.random()*mu_range*2,#-mu_range,
                           10,#random.random()*sigma_range
                           G)
             for i in verbose_gen(xrange(trials))]
    normalized_psums = [(pred,obs/pred) for pred,obs in psums]
    plt.scatter(*transpose(normalized_psums))
    m = max(concat(map(list,psums)))
    plt.plot([0,m],[0,m])
    plt.xlabel("Predicted LogZ")
    plt.ylabel("Experimental LogZ")
    return normalized_psums

def metaplot_powersum_test():
    """
    Does approximation converge as G -> \infty?
    """
    all_psums = []
    for i in range(1,4+1):
        print "starting on",i
        plt.subplot(2,2,i)
        all_psums.append(plot_powersum_test(G=1000*int(10**(i-1)),trials=10000/int(10**(i-1))))
    plt.show()
    return all_psums

def sum_lognormal(mu,sigma,n):
    xs = [random.gauss(mu,sigma) for _ in xrange(n)]
    y = sum([exp(x) for x in xs])
    

print "loaded"
