# The Makefile is one of the most important parts of the reproducible
# research workflow.

# Each component of the project (e.g. the final pub pdf, a png
# appearing in the pub, a csv containing data discussed in the paper,
# &c.) should be described in the makefile.  The makefile format is as
# follows:

# target_file: dependency1 (dependency2 dependency3 &c.)
# 	commands to generate target_file

# It is necessary to use tabs to indent.

results/coevolved_motifs/estremo_motifs/estremo_motifs.txt: data/estremo_motifs/alpha_sweep_motifs.py src/generate_coevolved_motifs/parse_estremo_motifs.py
	python src/generate_coevolved_motifs/parse_estremo_motifs.py results/coevolved_motifs/estremo_motifs/estremo_motifs.txt
