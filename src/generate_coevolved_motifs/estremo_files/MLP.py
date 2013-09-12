from Recognizer import Recognizer
from estremo_utils import transpose,beta,sigmoid,dot,takeBy,subst,sorted_indices,choose2,mmap,zipWith
from math import sqrt

class MLP(Recognizer):
    """Emulate an MLP"""
    
    @classmethod
    def validate_input(cls,recognizer_data):
        # drop the last line, which contains weights for the hidden layer
        input_weights = recognizer_data[:-1]
        hidden_weights = recognizer_data[-1]
        num_hidden_weights = len(hidden_weights)
        correct_input_line_length = all([len(line) % (4*num_hidden_weights) == 0
                                         for line in input_weights])
        correct_input_line_number = (len(input_weights) == num_hidden_weights)
        return correct_input_line_length and correct_input_line_number

    def __init__(self,recognizer_data):
        """Initialize an SLP"""
        # recognizer_data is of form:
        # [[XA1_weight,XT1_weight,XC1_weight,XG1_weight,XA2_weight,XT2_weight,XC2_weight,XG2_weight],
        #  [YA1_weight,YT1_weight,YC1_weight,YG1_weight,YA2_weight,YT2_weight,YC2_weight,YG2_weight],
        # ...
        #  [X_weight, Y_weight, Z_weight,..]                                                          ]
        # (for positions 1,2,3,... and hidden nodes X, Y, Z...)
        self.input_weights = [takeBy(4,hidden_node_weight)
                              for hidden_node_weight in recognizer_data[:-1]] #split by position
        self.hidden_weights = recognizer_data[-1]
        self.length = len(self.input_weights[0])

    def activation(self,site):
        input_activations = [sigmoid(sum([input_weight[i][Recognizer.base_index[base]]
                                 for i,base in enumerate(site)]))
                            for input_weight in self.input_weights]
        #print input_activations
        output_activation = sigmoid(dot(input_activations,self.hidden_weights))
        return output_activation
        
    def get_consensus(self):
        return "".join(Recognizer.inv_base_index[max(range(3),
                                                     key=lambda ind:col[ind])]
                       for col in self.input_weights)
        
                                

    def recognizer_asymmetry(self):
        # Implementing this for compatibility w/ SLP
        return 0

    def directional_score(self,site):
        return self.activation(site)

    def detect_correlations(self,site):
        """Does the MLP exhibit any correlations of the type
        demonstrated by Man & Stormo (NAR 2001)?"""
        bases = "ATCG" # follow ordering defined in Recognizer class
        w = len(site)
        all_dependent = True
        dependencies = 0
        discrepancy_dict = {}
        for (i,j) in choose2(range(w)):
            print (i,j)
            base1 = site[i]
            base2 = site[j]
            base_index1 = Recognizer.base_index[base1]
            base_index2 = Recognizer.base_index[base2]
            default_energy = self.binding_energy(site)
            alt_sites = [[subst(subst(site,b1,i),b2,j)
                          for b2 in bases]
                         for b1 in bases]
            alt_energies = mmap(self.binding_energy,alt_sites)
            best_energy = min(map(min,alt_energies))
            norm_energies = mmap(lambda x:x/best_energy,alt_energies)
            norm_default_energy = default_energy/best_energy
            independent = all((sorted_indices(alt_energy_row) ==
                               sorted_indices(alt_energies[0]))
                              for alt_energy_row in alt_energies)
            
            #independent = False
            if independent:
                all_dependent = False
            else:
                dependencies += 1
                ###Print motif
                print (i,j)
                print "     "
                for c2 in bases:
                    print "    ",c2,
                print
                
                for b_i1,c1 in enumerate(bases):
                    print c1,
                    for b_i2 in range(4):
                        print "%1.4f" % norm_energies[b_i1][b_i2],
                    print
                ###Finish printing motif
                base_1_mutant_energies = [(b,row[base_index2])
                                          for (b,row) in zip(bases,norm_energies)]
                print "base 1 mutants:",base_1_mutant_energies
                best_responses = [min(zip(bases,row),
                                      key=lambda (b,x):abs(norm_default_energy - x))
                                       for row in norm_energies]
                print "best responses:",best_responses
                discrepancy = max(zip(base_1_mutant_energies, best_responses),
                                  key = lambda ((b1,x),(b2,y)): abs(x-y))
                print "actual bases:",base1,base2,(base_index1,base_index2)
                print "norm default energy:",norm_default_energy
                (b1,v1),(b2,v2) = discrepancy
                print "discrepancy:",b1,b2,v1,v2,v2-v1
                discrepancy_dict[(i,j)] = v2-v1
                
        print "all dependent:",all_dependent
        print "# dependencies:", dependencies, "out of:",len(choose2(range(w)))
        return discrepancy_dict

    def detect_correlations2(self,site):
        """Does the MLP exhibit any correlations of the type
        demonstrated by Man & Stormo (NAR 2001)?"""
        bases = "ATCG" # follow ordering defined in Recognizer class
        w = len(site)
        all_dependent = True
        dependencies = 0
        discrepancy_dict = {}
        for (i,j) in choose2(range(w)):
            print (i,j)
            base1 = site[i]
            base2 = site[j]
            base_index1 = Recognizer.base_index[base1]
            base_index2 = Recognizer.base_index[base2]
            default_energy = self.binding_energy(site)
            alt_sites = [[subst(subst(site,b1,i),b2,j)
                          for b2 in bases]
                         for b1 in bases]
            alt_energies = mmap(self.binding_energy,alt_sites)
            best_energy = min(map(min,alt_energies))
            norm_energies = mmap(lambda x:x/best_energy,alt_energies)
            norm_default_energy = default_energy/best_energy
            independent = all((sorted_indices(alt_energy_row) ==
                               sorted_indices(alt_energies[0]))
                              for alt_energy_row in alt_energies)
            
            #independent = False
            if independent:
                all_dependent = False
            else:
                dependencies += 1
                ###Print motif
                print (i,j)
                print "     "
                for c2 in bases:
                    print "    ",c2,
                print
                
                for b_i1,c1 in enumerate(bases):
                    print c1,
                    for b_i2 in range(4):
                        print "%1.4f" % norm_energies[b_i1][b_i2],
                    print
                print "actual bases:",base1,base2,(base_index1,base_index2)
                ###Finish printing motif
                for base,row in zip(bases,norm_energies):
                    if (norm_energies[base_index1][base_index2] > row[base_index2] and
                        (all(norm_energies[base_index1][b] < row[b]
                             for b in range(4) if b != base_index1))):
                        print "discrepancy:", (norm_energies[base_index1][base_index2] -
                                              row[base_index2])
                        print "consider row:",base,row
                tr_norm_energies = transpose(norm_energies)
                for base,col in zip(bases,tr_norm_energies):
                    if (norm_energies[base_index1][base_index2] > col[base_index1]
                        and
                        (all(norm_energies[b][base_index2] < col[b]
                             for b in range(4) if b != base_index1))):
                        print "discrepancy:",(norm_energies[base_index1][base_index2] -
                                              col[base_index1])
                        print "consider col:",base,col
                      
                
        print "all dependent:",all_dependent
        print "# dependencies:", dependencies, "out of:",len(choose2(range(w)))
        return discrepancy_dict

    
    def l2_norm(self,other):
        """Compute L2 norm of the weights"""
        input_layer_sum = sum([(x - y)**2
                               for (self_hidden_node,
                                    other_hidden_node) in zip(self.input_weights,
                                                             other.input_weights)
                               for self_col,other_col in zip(self_hidden_node,
                                                             other_hidden_node)
                               for (x,y) in zip(self_col,other_col)])
        hidden_layer_sum = sum([(x - y)**2 for x,y in zip(self.hidden_weights,
                                                           other.hidden_weights)])
        return sqrt(input_layer_sum + hidden_layer_sum)
