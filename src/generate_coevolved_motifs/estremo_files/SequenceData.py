from estremo_utils import wc

class SequenceData(object):
    """Represent a site-describing line of ESTReMo output"""
    def __init__(self,seq_line):
        fields = seq_line.split(",")
        numeric_fields = map(float,fields[1:])
        self.site = fields[0]
        (self.position,
         self.strand,
         self.ep,
         self.act,
         self.tProb,
         self.gc,
         self.tConc) = numeric_fields[:7]
        self.position_energies = map(float,fields[6:])
        self.complements = sum(map(lambda(x,y):x==y,zip(self.site,wc(self.site))))

    def matched_base_pairs(self):
        return sum(map(lambda(x,y):x==y,
                       zip(self.site,wc(self.site))))

