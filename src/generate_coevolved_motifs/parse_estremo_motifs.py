import sys,os
print os.getcwd()
print sys.path
sys.path.append("data/estremo_motifs")
sys.path.append("src/generate_coevolved_motifs/estremo_files")
from Organism import Organism
print sys.path
from alpha_sweep_motifs import texts
from estremo_utils import random_site

orgs = map(lambda text:Organism(text,random_site(100000),overshooting=False),
               texts)
motifs = [org.motif for org in orgs]

if __name__ == "__main__":
    outfile = sys.argv[1]
    motif_strings = map("\n".join,motifs)
    with open(outfile,'w') as f:
        f.write("\n\n".join(motif_strings))
