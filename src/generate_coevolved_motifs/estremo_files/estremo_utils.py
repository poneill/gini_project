from __future__ import division
from math import log, exp,factorial,sqrt
import string, itertools, re, random
import scipy.misc as sm
from collections import Counter

Delta = "ATGC"
epsilon = 1e-10

kB  = 0.0019872041 #kcal/mol (!)
celsius_temp = 30
temp = celsius_temp + 273.15
beta = 1/(kB*temp) #NB: used to be defined as negative
estremo_beta = 1.0 / ((celsius_temp + 273.15) * 0.001987204118)

def verbose_gen(xs,modulus = 1):
    for i,x in enumerate(xs):
        if i % modulus == 0:
            print i
        yield x

def zipWith(f,xs,ys):
    return [f(x,y) for (x,y) in zip(xs,ys)]

def mean(xs):
    return sum(xs)/len(xs)

def median(xs):
    n = len(xs)
    if n % 2 == 1:
        return xs[n/2]
    else:
        return mean(xs[n/2-1:n/2+1])
    
def variance(xs):
    return mean(map(lambda x:x**2,xs)) - mean(xs)**2

def mean_and_stderr(xs):
    mu = mean(xs)
    stderr = sqrt(variance(xs)/len(xs))
    return mu,stderr

def log2(x):
    return log(x,2) if x > 0 else 0

def h(seq):
    return -1*sum([freq(b,seq)*log(freq(b,seq),2) for b in set(seq)]) 

complement_dict = {"a":"t","t":"a","g":"c","c":"g",
                   "A":"T","T":"A","G":"C","C":"G"}
def wc(seq):
    return "".join(complement_dict[c] for c in seq[::-1])

def freq(b,seq):
    return len([c for c in seq if c ==b])/len(seq)

def all_seqs_of_len(n,alphabet):
    if n == 1:
        return [a for a in alphabet]
    else:
        almost_all_seqs = all_seqs_of_len(n-1,alphabet)
    return flatten([map(lambda x: a+x,almost_all_seqs) for a in alphabet])

def equiprobable_seq_of_len(n,alphabet):
    return "".join(random.choice(alphabet) for i in range(n))

def partial_sums(alphabet_with_freqs):
    return [sum([awf[1] for awf in alphabet_with_freqs][:i+1]) for i in range(len(alphabet_with_freqs))]

def last(lst):
    return lst[-1]
def first(lst):
    return lst[0]

def random_seq_of_len(n,alphabet_with_freqs):
    """Return a sequence of length n with elements drawn from alphabet given as
    pairs of element, probability, e.g.: (("a",.25),("g",.5),...)"""
    elems = [awf[0] for awf in alphabet_with_freqs]
    cum_freqs = partial_sums(alphabet_with_freqs)
    elems_freqs = zip(elems, cum_freqs)
    def select(e_f):
        x = random.random()
        return first([elem for (elem, freq) in e_f if x < freq])
    return [select(elems_freqs) for i in range(n)]

def mod_out(seq):
    counts = sorted([seq.count(b) for b in 'atgc'])
    return [sum(col) for col in zip(counts, [0]+counts, [0,0]+counts)[:3]]

def product(ps):
    return reduce(lambda x,y:x*y,ps)

def group(xs):
    return sorted([[x for x in xs if x == s] for s in set(xs)])

def takeBy(n,xs):
    """Split xs into lists of size n"""
    return [xs[q:q+n] for q in range(0,len(xs),n)]
    
def mod_multiplicities(seq):
    n = len(seq)
    counts = sorted([seq.count(b) for b in 'atgc'])
    print counts
    mults = [len(g) for g in group(counts)]
    print mults
    nonempty = [count for count in counts if count > 0]
    print nonempty
    return (factorial(4)/product([factorial(m) for m in mults]) * (factorial(n) / 
            product([factorial(count) for count in counts])))

def mult_check(seq):
    return len(seqs_with_h_of(seq))

def seqs_with_h_of(seq):
    hseq = h(seq)
    epsilon = .000001
    return [s for s in all_seqs_of_len(len(seq), Delta) if abs(h(s) - hseq) < epsilon]

def check_ratio(seq):
    return mod_multiplicities(seq)/mult_check(seq)

def random_awf(elems):
    len_elems = len(elems)
    raw_freqs = sorted(random.random() for i in range(len_elems - 1))
    shift = zip(raw_freqs + [1], [0] + raw_freqs)
    freqs = [s[0] - s[1] for s in shift]
    return zip(elems, freqs)

def disc_freqs(input_string,n):
    d = {}
    for seq in all_seqs_of_len(n,set(input_string)):
        d[seq] = 0
    for i in range(len(input_string)-n+1):
        substr = input_string[i:i+n]
        d[substr] += 1
    return d

def disp_freqs(d):
    total = sum(d[key] for key in d)
    for key in sorted(d.keys()):
        print key, d[key]/total
        
def permutations(seq):
    if len(seq) == 1:
        return [seq]
    else:
        return flatten([map(lambda p:s+p,permutations(remove_first(s,seq))) for s in seq])

def remove_first(elem,lst):
    i = lst.index(elem)
    return lst[:i] + lst[i+1:]

def flatten(lst):
    return sum(lst,[])

def avg_h_of(seqs):
    return sum([h(seq) for seq in seqs])/len(seqs)

def htotal(k,n):
    return -(k/n*log2(k/n)+(n-k)/n*log2((n-k)/n)) * sm.comb(n,k)

def columnwise_h(seqs,correct=True):
    correction = 3/(2*log(2)*len(seqs)) if correct else 0 #Basharin 1959
    return [h(col) + correction for col in zip(*seqs)]

def mi(xs,ys):
    hx  = h(xs)
    hy  = h(ys)
    hxy = h(zip(xs,ys))
    return hx + hy - hxy

def columnwise_ic(seqs,hg=2,correct=True):
    return [hg - col_h for col_h in columnwise_h(seqs,correct=correct)]

def columnwise_mi(seqs, i, j):
    xs = [seq[i] for seq in seqs]
    ys = [seq[j] for seq in seqs]
    return mi(xs,ys)

def avg_h(n):
    return sum([htotal(k,n) for k in range(1,n)])/(2**n)

def r_freq(sites, genome):
    return log(len(genome) - len(sites))

def r_seq(sites,genome):
    hg = h(genome)
    return sum([hg - col_h for col_h in columnwise_h(sites)])

def error(n,genome):
    hg = h(genome)
    
# with open("t7genome.txt") as f:
#     genome = "".join(f.readlines()).strip()
    
# y= "[ct]"
# r= "[ag]"
# recognition_site = r"gt"+y+r+"ac"
# sites = re.findall(recognition_site,genome)
# hg = h(genome)

def instance_of_pattern(pattern):
    groups = ["".join(group) for 
              group in re.findall('\[(.*?)\]|([atgc])',pattern)]
    return "".join(random.choice(group) for group in groups)

def random_motif(n, length, alphabet_with_freqs):
    print "random motif"
    return [random_seq_of_len(length, alphabet_with_freqs) for i in range(n)]

def equiprobable_motif(n, length):
    """returns a list of n random sequences of specified length"""
    return [equiprobable_seq_of_len(length,Delta) for i in range(n)]

def non_random_motif(n, pattern):
    return [instance_of_pattern(pattern) for i in range(n)]

#motifs = [non_random_motif(100,"at[at][gc]gc") for i in range(100)]

def gini(xs):
    ys = sorted(xs)
    n = float(len(ys))
    return (2*sum((i+1)*y for i,y in enumerate(ys)))/(n*sum(ys)) - (n+1)/n

def motif_gini(motif):
    """Return the gini coefficient of the column ics"""
    return gini(columnwise_ic(motif))
    
    
# col_icx = []
# col_icy = []

def get_ics_mis_from_motifs(motifs):
    print "made it this far"
    k = 0
    avg_col_ics = []
    pairwise_mis = []
    for motif in motifs:
        print k
        k+=1
        colwise_ics = columnwise_ic(motif)
        for i in range(len(motif[0])):
            for j in range(i+1,len(motif[0])):
                # col_icx.append(colwise_ics[i])
                # col_icy.append(colwise_ics[j])
                avg_col_ic = (colwise_ics[i] + colwise_ics[j])/2
                pairwise_mi = columnwise_mi(motif,i,j)
                avg_col_ics.append(avg_col_ic)
                pairwise_mis.append(pairwise_mi)
    return avg_col_ics, pairwise_mis

def string_replace(st,index,c):
    return st[:index] + c + st[index + 1:]

def motif_ic(motif,correct=True):
    return sum(columnwise_ic(motif,correct=True))

def random_site(n):
    return "".join([random.choice("ACGT") for i in range(n)])

def random_palindrome(length):
    if length % 2 == 0:
        half_site = random_site(length//2)
        return half_site + wc(half_site)
    else:
        half_site = random_site((length-1)//2)
        middle = random_site(1)
        return half_site + middle + wc(half_site)

def random_palindromic_motif(n,length):
    return [random_palindrome(length) for i in range(n)]

def generate_greedy_motif_with_ic(desired_ic,epsilon,num_seqs,length,verbose=False):
    motif = [random_site(length) for i in range(num_seqs)]
    ic = motif_ic(motif)
    while(abs(desired_ic - ic) > epsilon):
        motif_prime = motif[:]
        n = random.randrange(num_seqs)
        l = random.randrange(length)
        motif_prime[n] = string_replace(motif_prime[n],l,random.choice("ACGT"))
        ic_prime = motif_ic(motif_prime)
        if abs(ic_prime - desired_ic) < abs(ic - desired_ic):
            motif = motif_prime
            ic = motif_ic(motif)
        if verbose:
            print ic
    return motif

def control_motif(motif,epsilon=0.01,verbose=False):
    return generate_greedy_motif_with_ic(motif_ic(motif),epsilon,
                                         len(motif),len(motif[0]),verbose)
    
    
print("loaded")
    
    
def subsample_motif(motif,n):
    select = [random.randrange(len(motif)) for i in range(n)]
    return [motif[s] for s in select]

def crp_vs_control_test(trials,motif_size,filename):
    data = []
    for i in range(trials):
	crp_boot = subsample_motif(crp,motif_size)
	ic = motif_ic(crp_boot)
 	control = generate_greedy_motif_with_ic(ic,.1,motif_size,22)
 	g,g_control = gini(columnwise_ic(crp_boot)),gini(columnwise_ic(control))
        print g, g_control
        data.append([g,g_control])
    data2csv(data,filename)

def transpose(xxs):
    return map(list,zip(*xxs))

def col_prob(col):
    return [col.count(c)/float(len(col)) for c in Delta]

def kl(ps,qs):
    return sum(p*log2(p/(q+epsilon)) for (p,q) in zip(ps,qs))

def motif_kl(motif1,motif2):
    assert(len(motif1[0]) == len(motif2[0]))
    cols1,cols2 = transpose(motif1),transpose(motif2)
    total_kl = sum(kl(col_prob(col1),col_prob(col2))
                   for (col1,col2) in zip(cols1,cols2))
    return total_kl

def wckl(motif):
    return motif_kl(motif,map(wc,motif))

def sigmoid(x):
    return 1/(1 + exp(-x))

def dot_reference(xs,ys):
    """This implementation is apparently too cool for python.  see dot"""
    return sum(zipWith(lambda x,y:x*y,xs,ys))

def dot(xs,ys):
    """Fast dot product operation"""
    acc = 0
    for (x,y) in zip(xs,ys):
        acc += x * y
    return acc

def invert_dict(d):
    return {v:k for (k,v) in d.items()}

def subst(xs,c,i):
    """Return xs with c replacing xs[i]"""
    if type(xs) is str:
        return xs[:i] + c + xs[i+1:]
    else:
        return xs[:i] + [c] + xs[i+1:]

def choose2(xs,gen=False):
    """return list of choose(xs, 2) pairs, retaining ordering on xs"""
    if gen:
        return ((x1, x2) for i, x1 in enumerate(xs) for x2 in xs[i+1:])
    else:
        return [(x1, x2) for i, x1 in enumerate(xs) for x2 in xs[i+1:]]

def mmap(f,xxs):
    return [map(f,xs) for xs in xxs]

def sorted_indices(xs):
    """Return a list of indices that puts xs in sorted order.
    E.G.: sorted_indices([40,10,30,20]) => [1,3,2,0]"""
    return [i for (i,v) in sorted(enumerate(xs),key=lambda(i,v):v)]

def test_sorted_indices(xs):
    si = sorted_indices(xs)
    return sorted(xs) == [xs[i] for i in si]

def get_ecoli_genome(at_lab=False):
    lab_file = "/home/poneill/ecoli/NC_000913.fna"
    home_file = "/home/pat/ecoli/NC_000913.fna"
    with open(lab_file if at_lab else home_file) as f:
        genome = "".join([line.strip() for line in f.readlines()[1:]])
    return "".join(g for g in genome if g in "ATGC") # contains other iupac symbols

def random_substring(xs,k):
    i = random.randint(0,len(xs)-k)
    return xs[i:i+k]

def nmers(genome,n):
    """Return a sliding window of width n"""
    for i in xrange(len(genome)- n + 1):
        yield genome[i:i+n]
    
def nmer_dict(genome,n,double_stranded=False):
    """Return a dictionary of counts for each nmer in genome"""
    fd_counts = Counter(nmers(genome,n))
    if not double_stranded:
        return fd_counts
    else:
        bk_counts = Counter(nmers(wc(genome),n))
        return fd_counts + bk_counts

def concat(xxs):
    return sum(xs,[])

def normalize(xs):
    z = float(sum(xs))
    return [x/z for x in xs]
