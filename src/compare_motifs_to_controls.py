"""
This script generates boxplots of motif-recognizer vs. control motif
Gini coefficients.

Inputs: a set of motifs and a selection variable for the desired type of control

Outputs: a boxplot comparing the Gini coefficient of the observed motifs to those of an ensemble of matched random controls.

"""

def make_gini_comparison_plot(motifs,control_ensembles,descriptor,filename=None):
    motif_ginis = map(motif_gini,motifs)
    indices = sorted_indices(motif_ginis)
    sorted_motif_ginis = rslice(motif_ginis,indices)
    sorted_control_ensembles = rslice(control_ensembles,indices)
    plt.boxplot(sorted_control_ensembles,label="Control Ensemble")
    plt.scatter(range(1,len(motifs) + 1),sorted_motif_ginis,
             label="Observed")
    plt.ylabel("Gini Coefficient")
    plt.legend(loc=2)
    title = "Gini Coefficients for %s Motifs and Random Ensembles" % descriptor
    plt.title(title)
