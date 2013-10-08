from motifs import *
from utils import *
from matplotlib import pyplot as plt

ecoli_motifs = [getattr(Escherichia_coli,tf)
               for tf in Escherichia_coli.tfs]
fly_motifs = [getattr(Drosophila_melanogaster,tf)
               for tf in Drosophila_melanogaster.tfs]

ecoli_ginis = map(motif_gini,ecoli_motifs)
ecoli_ics = map(motif_ic,ecoli_motifs)

fly_ginis = map(motif_gini,fly_motifs)
fly_ics = map(motif_ic,fly_motifs)

plt.scatter(ecoli_ics,ecoli_ginis)
plt.scatter(fly_ics,fly_ginis,color='g')
plt.xlabel("IC")
plt.ylabel("Gini coefficient")
plt.show()
