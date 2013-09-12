from math import log,exp
from estremo_utils import beta,wc,invert_dict
import random

class Recognizer(object):
    """Wrapper class for SLP, MLP"""
    base_index = {
        "A":0,
        "T":1,
        "C":2,
        "G":3
        }
    inv_base_index = invert_dict(base_index)
    
    def __init__(self,recognizer_data):
        """Should not be called: use the RecognizerFactory instead.
        We are doing it this way in order to avoid circular dependency
        issues."""
        pass

    def binding_energy(self,site,double_stranded=False):
        """Compute binding energy of site.  See calcBindingEnergy,Organism.cpp"""
        w = len(site)
        E_f = self.activation(site) * w * (-2/beta) - 0.05
        E_b = double_stranded and (self.activation(wc(site)) * w * (-2/beta) - 0.05)
        E_t = ((1 / -beta) * log(exp(-beta * E_f) + exp(-beta * E_b))
               if double_stranded
               else E_f)
        return E_t

    def score(self,site,both_strands=False):
        ef = self.directional_score(site)
        if both_strands:
            eb = self.directional_score(wc(site))
            return log(exp(beta*ef) + exp(beta*eb))/(beta)
        else:
            return ef

    def binding_energy_statistic(self,f=min,num_trials=10000):
        """Return some statistic of ensemble of random binding energies"""
        return f([self.binding_energy(site,double_stranded=True)
                    for site in ["".join([random.choice("ATGC")
                                          for i in range(self.length)])
                                 for i in xrange(num_trials)]])

    
