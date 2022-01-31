## mycomimic
# Randomly generate synthetic ITS2 sequences for safe and easy use in fungal metabarcoding projects

While synthetic DNA spike-ins derived from gene fragments can help researchers obtain more quantitative data from their metabarcoding efforts,
doing requires cumbersome optimzation of the amount of input spike-in to add to each sample and still requires sacrificing
a substantial portion of sequence reads to obtain a sufficiently reliable quantitative adjustment.

The MycoMimic approach, while still in development, is designed to address both of these concerns. Users don't have to
get the input amount of spike-in just right for their samples and instead can use blocking oligonucleotides to reduce
the spike-in amplification at PCR in the event that the input amount is overshot. Spike-in and natural fungal amplicons
can then be visualized separately for each sample after a restriction enzyme digest.

While the MycoMimic method still requires read sacrifice, these reads can be used to identify and correct for the bias
among degenerate primers occurring in each sample. By ligating oligonucleotide duplexes to the main ITS2-containing gene
fragment, the priming regions of each spike-in can have degenerate bases. Researchers can then bioinformatically tally up
the sequence counts of each primer variant in each sample to correct for the over- or under-abundance of each variant,
potentially using these tallies to inform a correction model for every sample in the sequencing effort.
