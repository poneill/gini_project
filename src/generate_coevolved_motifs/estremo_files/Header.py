class Header(object):
    def __init__(self,header_line):
        fields = [field.strip() for field in header_line.split(",")]
        (self.gen,
         self.r_c,
         self.r_f,
         self.r_s,
         self.mi,
         self.gGC,
         self.nGC,
         self.sumBG,
         self.sumNCR,
         self.numTF,
         self.pb_err,
         self.chemPot,
         self.ncr_mu,
         self.rec_mu,#recognizer mutation
         self.beta,
         self.bg,
         self.fitness,
         self.maxFitness)= map(float,fields[:-1]) #all but last field
        self.parent = int(fields[-1],16) # last field is in hex; convert it to dec.
        
