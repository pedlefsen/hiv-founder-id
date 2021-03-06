# hiv-founder-id
Estimate founders of HIV-1 infection using nucleotide sequences from very acute infection.  This consists of a handful of perl scripts with one main entry point: identify-founders.pl.  These scripts use the InSites program (via Jim Mullins's lab website, DIVEIN, at http://indra.mullins.microbiol.washington.edu/DIVEIN/insites.html) to process an input nucleotide multiple sequence alignment of HIV-1 sequences from a single host at an early timepont to calculate so-called informative and private sites: columns of the alignment in which at least two residues (or gap) each appear in at least two sequences are called "informative", and columns in which exactly one sequence differs from the consensus are called "private".  The software computes the ratio of informative to private sites (higher values indicate greater departure from what is expected under a single-founder Poisson model), and also clusters the sequences and computes consensus sequences for each cluster.  These are the estimated founders of the HIV-1 infection (when the informative-to-private site ratio is hgih).  This procedure does not evaluate recombination or hypermutation at this time.

Dependencies:

ln -s <path to patched phyml executable> phyml
ln -s <path to profuse> profuse
