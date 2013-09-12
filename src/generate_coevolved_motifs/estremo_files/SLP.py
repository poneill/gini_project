from estremo_utils import transpose, wc, beta
from math import exp,log,sqrt
from Recognizer import Recognizer

class SLP(Recognizer):
    """Emulate an SLP"""

    @classmethod
    def validate_input(cls,recognizer_data):
        """Are there four lines, one corresponding to each base,
        containing weights the input nodes?"""
        return len(recognizer_data) == 4

    def __init__(self,recognizer_data):
        """Initialize an SLP."""
        #input layer is of form [[A1_weight,A2_weight,A3_weight,...,Aw_weight],
        #                        [T1_weight,T2_weight,T3_weight,...,Tw_weight],
        #                        [C1_weight,...                              ],
        #                        [G1_weight,...                              ]]

        self.input_layer = recognizer_data
        self.length = len(recognizer_data[0])
        
    def recognizer_asymmetry(self):
        A = 0
        T = 1
        C = 2
        G = 3
        positions = transpose(self.input_layer)
        end = len(positions) - 1
        asymmetry = sum([((positions[i][A] - positions[end-i][T])**2 +
                          (positions[i][T] - positions[end-i][A])**2 +
                          (positions[i][G] - positions[end-i][C])**2 +
                          (positions[i][C] - positions[end-i][G])**2)
                         for i in range(int((end + 1)/2))])
        return asymmetry

    def activation(self,site):
        positions = transpose(self.input_layer)
        total = sum([positions[i][Recognizer.base_index[base]]
                     for i,base in enumerate(site)])
        return 1/(1+exp(-total))

    def directional_score(self,site):
        return self.activation(site)
    
    def motif_sites(self):
        """Return a motif that approximates, in some sense, the
        activity of the recognizer"""
        motif = []
        while (len(motif) < 100):
            site = "".join([random.choice("ATCG") for i in range(self.length)])
            s = self.score(site)
            if random.random() < s:
                motif.append(site)
        return motif

    def l2_norm(self,other):
        """Compute L2 norm of the weights"""
        return sqrt(sum([(x - y)**2
                             for self_col,other_col in zip(self.input_layer,
                                                           other.input_layer)
                             for (x,y) in zip(self_col,other_col)]))

    def min_binding_site(self):
        """Return the site yielding the minimum binding energy"""
        return "".join([min(zip("ATCG",col),key=lambda (c,w):w)[0]
                        for col in transpose(self.input_layer)])

    def max_binding_site(self):
        """Return the site yielding the minimum binding energy"""
        return "".join([max(zip("ATCG",col),key=lambda (c,w):w)[0]
                        for col in transpose(self.input_layer)])
