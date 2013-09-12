from Paragraph import Paragraph
from estremo_utils import beta,random_site,string_replace,random_substring,get_ecoli_genome,mean,nmer_dict
from math import exp

# check overshooting

class Organism(object):
    """Emulate an ESTReMo-style organism"""
    def __init__(self,paragraph_text,background_string,overshooting=True,
                 background_scaling_factor=100,
                 ncr_size=32):
        """Initialize an organism from a paragraph and a background"""
        paragraph_lines = paragraph_text.split("\n")
        self.paragraph = Paragraph(paragraph_lines)
        self.targets = [seq.tProb for seq in self.paragraph.seq_data]
        self.motif = self.paragraph.motif()
        self.original_motif = self.motif[:] #copy for backup against mutations
        self.recognizer = self.paragraph.recognizer
        self.background = background_string
        self.memo_background_Z = None
        self.overshooting = overshooting
        self.background_scaling_factor = background_scaling_factor
        self.ncr_size = 32
        # samples * scaling factor = effective background
        
    def fitness(self):
        site_energies = [self.recognizer.binding_energy(site)
                         for site in self.motif]
        w = len(self.motif[0])
        n = len(self.background)
        background_energies = [self.recognizer.binding_energy(self.background[i:i+w])
                                                   for i in range(n-w+1)]
        foreground_Z = sum([exp(-beta * energy) for energy in site_energies])
        Z = foreground_Z + self.background_Z()
        site_occupancies = [exp(-beta*energy)/Z for energy in site_energies]
        print "site occupancies:",site_occupancies
        return sum(self.site_fitness(occ,target)
                   for (occ,target) in zip(site_occupancies,self.targets))

    def background_Z(self):
        if not self.memo_background_Z is None:
            return self.memo_background_Z
        else:
            w = len(self.motif[0])
            n = len(self.background)
            samples = n / self.background_scaling_factor
            print "samples:",samples
            be = self.recognizer.binding_energy
            bkgd_energies = [be(random_substring(self.background,w))
                             for i in range(samples)]
            self.memo_background_Z = sum([exp(-beta * energy)
                                          for energy in bkgd_energies]) * self.background_scaling_factor
        return self.memo_background_Z
    
    def site_fitness(self,occupancy,target):
        # print "occupancy:",occupancy
        # print "target:",target
        delta = 0.17
        eta = 0.02
        M = 1.8
        Ky = 0.4  #Params come from config file
        Z = occupancy * M
        L = Ky / ((delta / eta) * (1 - target)**2 - 1)
        g = delta * ((Z * L) / (Ky + L)) - eta * (Z / (1 - (Z/M)))
        if not self.overshooting: #i.e. if overshooting is not penalized...
            # compute optimum expression level
            Zopt = M * (1 - sqrt((eta / delta) * ((L + Ky) / L)))
            # and fitness corresponding to optimum expression
            gopt = delta * ((Zopt * L) / (Ky + L)) - eta * (Zopt / (1 - (Zopt/M)))
            if Z > Zopt:
                return gopt
        #print g
        return g
    
    def mutate_site(self,site_number,position,base):
        site = self.motif[site_number]
        self.motif[site_number] = string_replace(site,position,base)

    def reset_motif(self):
        self.motif = self.original_motif[:]

    def explore_site_mutations(self):
        print self.motif
        self.original_motif = self.motif[:]
        fit = self.fitness()
        mutations = []
        for site_number in range(len(self.motif)):
            for position in range(len(self.motif[0])):
                original_base = self.motif[site_number][position]
                for base in "ACTG":
                    if original_base == base:
                        continue #don't recompute original fitness'
                    self.mutate_site(site_number,position,base)
                    fit_prime = self.fitness()
                    report_string = "Improvement" if fit_prime > fit else ""
                    diff = (fit_prime - fit)/fit
                    mutations.append((site_number,position,base,fit_prime))
                    print site_number,position,base,self.fitness(),diff,report_string
                self.reset_motif()
        return mutations

    def grad_descent(self):
        iteration = 0
        mutations = None
        motif_dustbin = []
        while (mutations is None or
               any(fit > 0 for (site,pos,base,fit) in mutations)):
            mutations = self.explore_site_mutations()
            site_idx, pos, base, fit = max(mutations,
                                           key=lambda tup:tup[3])#max by fitness
            self.motif[site_idx] = string_replace(self.motif[site_idx],pos,base)
            print "iteration:",iteration,site_idx,pos,base,fit
            motif_dustbin.append(self.motif[:])
            iteration += 1

    def serialize(self,population_size):
        """Print self out per ESTReMo's population serialization format"""
        print population_size
        for i in range(population_size):
            print i
            print len(self.motif)
            for site in self.motif:
                # TODO we should probably define the organism to
                # include the ncr.  For now, randomize the rest of the
                # ncr
                print site + random_site(self.ncr_size - len(site))
            print 1 if type(self.recognizer) is SLP else 2
            # TODO needs to be generalized to MLPs!
            weights = concat(map(list,transpose(self.recognizer.input_layer)))
            print len(weights)
            for weight in weights:
                print weight

print "loaded Organism"
