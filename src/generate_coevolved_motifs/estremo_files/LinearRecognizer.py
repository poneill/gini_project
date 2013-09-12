from Recognizer import Recognizer
from estremo_utils import transpose,beta,sigmoid,dot,takeBy,subst,sorted_indices,choose2,mmap,zipWith
from math import sqrt

class LinearRecognizer(Recognizer):
    """Emulate an MLP"""
    
    @classmethod
    def validate_input(cls,recognizer_data):
        # We forgot to add in newlines, so the linear recognizer has
        # its weights arrayed in one long line...
        return len(recognizer_data) == 4

    def __init__(self,recognizer_data):
        """Initialize a linear recognizer"""
        # recognizer_data is of form:
        # [[XA1_weight,XT1_weight,XC1_weight,XG1_weight,XA2_weight,XT2_weight,XC2_weight,XG2_weight],
        #  [YA1_weight,YT1_weight,YC1_weight,YG1_weight,YA2_weight,YT2_weight,YC2_weight,YG2_weight],
        # ...
        #  [X_weight, Y_weight, Z_weight,..]                                                          ]
        # (for positions 1,2,3,... and hidden nodes X, Y, Z...)
        #print recognizer_data
        self.length = len(recognizer_data[0])
        #print "length:",self.length
        self.input_layer = recognizer_data

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
        return total

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
