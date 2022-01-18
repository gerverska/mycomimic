# SYNTHETIC ITS2 GENE FRAGMENT GENERATION ####
# This gene fragment is compatible with the following fungal primers:
# fITS7, ITS3, ITS3_KYO1, 5.8S-Fun, ITS4, ITS4-Fun, ITS4ngs
# Take care to size ITS2 appropriately to your needs!

# This approach improves on J. Palmer's "amptk_synthetic_mock.py" by:
# > increasing the length of the spike-in to facilitate protocol development and downstream library re-balancing
# > generating sequences with accurate GC%
# > reducing the potential for sequence shuffling functions to put bases back where they started
# > adding blocking primer regions (annealing inhibition- and elongation arrest-type) to allow for post-extraction tweaking of spike-in presence
# > including the option to insert adapter sequences between universal primers to allow for the resizing of synthetic amplicons at the second stage of PCR

# Load packages ####
library(tidyverse)
library(magrittr)
library(here)
library(cowplot)

# Load functions ####
random.seq.gc <- function(bp, gc){
  
  # Determine the number of bases for each base group
  bp.gc <- bp %>% multiply_by(gc) %>% round()
  bp.at <- bp - bp.gc
  
  # For odd values of bp.at, bp.a < bp.t
  bp.a <- bp.at %>% divide_by(2) %>% floor()
  bp.t <- bp.at - bp.a
  
  # For odd values of bp.gc, bp.g > bp.c
  bp.g <- bp.gc %>% divide_by(2) %>% ceiling()
  bp.c <- bp.gc - bp.g
  
  # The sample() provides an initial shuffle
  c(replicate(bp.a, 'A'), replicate(bp.t, 'T'),
    replicate(bp.g, 'G'), replicate(bp.c, 'C')) %>% sample() %>% paste(collapse = '')
  
}
seq.shuffle <- function(string, swap.rate = 0.5, shuffles = 10){
  # Maximum swap.rate is 0.5--half of the bases swap with the other half
  # In cases where there is an odd number of bases, a random base is left unswapped.
  
  vect <- strsplit(string, '')[[1]]
  
  count <- 0
  
  while(shuffles - 1 >= 0){
    
    count <- count + 1
    
    if(swap.rate == 'random'){
      swap.rate <- runif(1, max = 0.5)
    }
    
    swaps <- vect %>% length() %>% multiply_by(swap.rate) %>% floor()
    swappers <- sample(1:length(vect), swaps)
    swappees <- setdiff(1:length(vect), swappers) %>% sample(swaps)
    
    swapper.bases <- vect[swappers]
    swappee.bases <- vect[swappees]
    
    vect[swappers] <- swappee.bases
    vect[swappees] <- swapper.bases
    
    shuffles <- shuffles -1
    
  }
  
  vect %>% paste(collapse = '')
  
}
gene.generator <- function(perm = 100, adapter = 'no', frameshift = 'no'){
  
  # Define ribosomal subunit DNA sequences, derived from Saccharomyces cerevisiae S288C
  # Both 5.8S and LSU are trimmed to only contain forward and reverse primers of interest
  upstream_5.8s <- 'AACTTTCAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGC' # 49 bp, GC content = 0.469
  
  fits7 <- 'GTGAATCATCGAATCTTTG' # 19 bp, GC content = 0.368
  
  its4.all <- 'ACTTAAGCATATCAATAAGCGGAGGA' # 26 bp, GC content = 0.385
  
  # Set sequence generation parameters
  upstream_5.8s.bp <- 49
  upstream_5.8s.gc <- 0.469
  mid_5.8s.bp <- 34
  mid_5.8s.gc <- 0.382
  fits7.bp <- 19
  fits7.gc <- 0.368
  downstream_5.8s.bp <- 55
  downstream_5.8s.gc <- 0.545
  its2.bp <- 500
  its2.gc <- 0.5
  upstream_28s.bp <- 34
  upstream_28s.gc <- 0.5
  its4.all.bp <- 26
  its4.all.gc <- 0.385
  
  # Request iterations
  perm <- 1:perm
  
  for(i in perm){
    
    # Generate random sequences with a set GC content
    # It's still not apparent to me how to optimize the number of shuffles...
    
    mid_5.8s <- random.seq.gc(mid_5.8s.bp, mid_5.8s.gc) %>%
      seq.shuffle(shuffles = 10)
    downstream_5.8s <- random.seq.gc(downstream_5.8s.bp, downstream_5.8s.gc) %>%
      seq.shuffle(shuffles = 10)
    
    if(adapter == 'forward' || adapter == 'reverse'){
      
      # If 'forward' or 'reverse' is provided as an argument, the below adapter sequences will be inserted
      # in the middle of a randomly generated ITS sequence, along with a specified number of frameshift Ns.
      # This forward or reverse options are provided to allow for situations in which the traditionally
      # 'forward' adapter is added to the traditionally 'reverse' universal primer (see Taylor et al. 2016, 5.8S-Fun + ITS4-Fun).
      # Choose 'reverse' if your primers have the conventional arrangement (forward adapter with forward primer, etc.).
      # In both cases, the reverse complement sequence is specified to ensure proper integration into the '-' strand.
      f.adapter.rc <- 'CTGTCTCTTATACACATCTGACGCTGCCGACGA' # 33 bp, GC content = 0.515
      r.adapter.rc <- 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC' # 34 bp, GC content = 0.529
      
      its2.part.bp <- its2.bp / 2
      its2.part.gc <- its2.gc
      
      its2.pt.1 <- random.seq.gc(its2.part.bp, its2.part.gc) %>%
        seq.shuffle(shuffles = 10)
      
      fs.seq <- rep('N', frameshift) %>% paste(collapse = '')
      
      its2.pt.2 <- random.seq.gc(its2.part.bp, its2.part.gc) %>%
        seq.shuffle(shuffles = 10)
      
      if(adapter == 'forward'){
        
        adapter.seq <- f.adapter.rc
        
      } else {
        
        adapter.seq <- r.adapter.rc
        
      }
      
      its2 <- paste0(its2.pt.1, fs.seq, adapter.seq, its2.pt.2)
      
    } else {
      
      adapter.seq <- 'no'
      fs.seq <- 'no'
      its2.part.bp <- 'no'
      its2.part.gc <- 'no'
      its2.pt.1 <- 'no'
      its2.pt.2 <- 'no'
      
      its2 <- random.seq.gc(its2.bp, its2.gc) %>%
        seq.shuffle(shuffles = 20)
      
    }
    
    upstream_28s <- random.seq.gc(upstream_28s.bp, upstream_28s.gc) %>%
      seq.shuffle(shuffles = 10)
    
    # Combine strings to produce the full, synthetic ITS
    mimic.name <- paste0('MycoMimic.', i)
    mimic.seq <- paste0(upstream_5.8s, mid_5.8s, fits7, downstream_5.8s, its2, upstream_28s, its4.all)
    
    # Prepare FASTA output for the full gene block, as well as separate outputs for each randomly generated component
    mimic.out <- paste0('>', mimic.name, ';adapter=', adapter, ';frameshift=', frameshift,
                        ';bp_ITS2=', its2.bp, ';GC_ITS2=', its2.gc, '\n',
                        mimic.seq)
    
    fasta.out <- paste0('mycomimic_', adapter, '.adapter.fasta')
    
    if(i == 1){
      
      mimic <- data.frame(gene = mimic.name,
                          upstream_5.8s.bp = upstream_5.8s.bp, upstream_5.8s.gc = upstream_5.8s.gc, upstream_5.8s = upstream_5.8s,
                          mid_5.8s.bp = mid_5.8s.bp, mid_5.8s.gc = mid_5.8s.gc, mid_5.8s = mid_5.8s,
                          fits7.bp = fits7.bp, fits7.gc = fits7.gc, fits7 = fits7,
                          downstream_5.8s.bp = downstream_5.8s.bp, downstream_5.8s.gc = downstream_5.8s.gc, downstream_5.8s = downstream_5.8s,
                          adapter = adapter, frameshift = fs.seq, adapter.seq = adapter.seq,
                          its2.pt.1.bp = its2.part.bp, its2.pt.1.gc = its2.part.gc, its2.pt.1 = its2.pt.1,
                          its2.pt.2.bp = its2.part.bp, its2.pt.2.gc = its2.part.gc, its2.pt.2 = its2.pt.2,
                          its2.bp = its2.bp, its2.gc = its2.gc, its2 = its2,
                          upstream_28s.bp = upstream_28s.bp, upstream_28s.gc = upstream_28s.gc, upstream_28s = upstream_28s,
                          its4.all.bp = its4.all.bp, its4.all.gc = its4.all.gc, its4.all = its4.all,
                          full = mimic.seq)
      
      
      
      cat(mimic.out, file = here(sequences, fasta.out), sep = '\n')
      
    } else {
      
      mimic.buffer <- data.frame(gene = mimic.name,
                                 upstream_5.8s.bp = upstream_5.8s.bp, upstream_5.8s.gc = upstream_5.8s.gc, upstream_5.8s = upstream_5.8s,
                                 mid_5.8s.bp = mid_5.8s.bp, mid_5.8s.gc = mid_5.8s.gc, mid_5.8s = mid_5.8s,
                                 fits7.bp = fits7.bp, fits7.gc = fits7.gc, fits7 = fits7,
                                 downstream_5.8s.bp = downstream_5.8s.bp, downstream_5.8s.gc = downstream_5.8s.gc, downstream_5.8s = downstream_5.8s,
                                 adapter = adapter, frameshift = fs.seq, adapter.seq = adapter.seq,
                                 its2.pt.1.bp = its2.part.bp, its2.pt.1.gc = its2.part.gc, its2.pt.1 = its2.pt.1,
                                 its2.pt.2.bp = its2.part.bp, its2.pt.2.gc = its2.part.gc, its2.pt.2 = its2.pt.2,
                                 its2.bp = its2.bp, its2.gc = its2.gc, its2 = its2,
                                 upstream_28s.bp = upstream_28s.bp, upstream_28s.gc = upstream_28s.gc, upstream_28s = upstream_28s,
                                 its4.all.bp = its4.all.bp, its4.all.gc = its4.all.gc, its4.all = its4.all,
                                 full = mimic.seq)
      mimic <- rbind(mimic, mimic.buffer)
      cat(mimic.out, file = here(sequences, fasta.out), sep = '\n', append = T)
      
    }
    
  }
  
  # Write a summary csv and output a reference data frame after calculating the total length of each gene block
  # mimic$length <- mimic$upstream_5.8s.bp + mimic$mid_5.8s.bp + mimic$fits7.bp + mimic$downstream_5.8s.bp +
  #   mimic$its2.bp +
  #   mimic$upstream_28s.bp + mimic$its4.all.bp
  write.csv(mimic, here(logs, 'mycomimic.csv'), row.names = F)
  
  return(mimic)
  
}
initiate.inhibitors <- function(df, primer){
  
  # Different parameters are associated with each primer.
  # Current starting indices correspond to the maximum overlap for each primer.
  # If the primer is on the minus strand, a reverse complement is generated.
  if(primer == 'fits7'){
    
    region <- df[, c('gene', 'fits7', 'downstream_5.8s')]
    region$sequences <- paste0(region[, 2], region[, 3])
    
    max.overlap <- 10
    min.overlap <- 1
    start <- 10:39
    rc <- F
    
  } else if(primer == 'its3.kyo1'){
    
    region <- df[, c('gene', 'upstream_5.8s', 'mid_5.8s')]
    region$sequences <- paste0(region[, 2], region[, 3])
    
    max.overlap <- 9
    min.overlap <- 1
    start <- 41:70 
    rc <- F
    
  } else if(primer == 'its3'){
    
    region <- df[, c('gene', 'upstream_5.8s', 'mid_5.8s')]
    region$sequences <- paste0(region[, 2], region[, 3])
    
    max.overlap <- 10
    min.overlap <- 1
    start <- 40:69
    rc <- F
    
  }
  
  else if(primer == 'its4.fun'){
    
    region <- df[, c('gene', 'upstream_28s', 'its4.all')]
    region$sequences <- paste0(region[, 2], region[, 3])
    
    max.overlap <- 13
    min.overlap <- 1
    start <- 18:47
    rc <- T
    
  }
  
  inhibitor.out.name <- paste0(primer, '_inhibitors.fasta')
  
  # Different 30-mer sequences (i.e., the blocking primer candidates) are generated here with decreasing amounts of overlap for the specified primer.
  # The number of available 30-mers depends on the length of the specified primer.
  # Candidates with greater amounts of overlap run a greater risk of blocking the amplification of real fungal amplicons,
  # while candidates with less overlap might not block enough synthetic amplification.
  # This approach allows us to maximize overlap and minimize--ideally, entirely avoid--similarity to real fungal priming regions.
  for(i in 1:nrow(df)){
    
    gene <- region[i, 'gene']
    vect <- strsplit(region[i, 'sequences'], '')[[1]]
    
    total.shifts <- max.overlap - min.overlap
    index.shifts <- 0:total.shifts
    
    for(j in index.shifts){
      
      if(rc == F){
        
        shifted.indices <- start + j
        shifted.seq <- vect[shifted.indices] %>% paste(collapse = '')
        
      } else {
        
        shifted.indices <- start - j
        shifted.seq <- vect[shifted.indices] %>% rev() %>% paste(collapse = '')
        shifted.seq <- chartr('ATGC','TACG', shifted.seq)
        
        shifted.indices <- shifted.indices %>% rev()
        
      }
      
      current.overlap <- max.overlap - j
      inhibitor.id <- paste(gene, primer, shifted.indices[[1]], current.overlap, sep = '_')
      inhibitor.out <- paste0('>', inhibitor.id, '\n',
                            shifted.seq)
      
      if(i == 1 && j == 0){
        
        inhibitor.vect <- inhibitor.id
        cat(inhibitor.out, file = here(inhibitors, inhibitor.out.name), sep = '\n')
        
      } else {
        
        inhibitor.buffer <- inhibitor.id
        inhibitor.vect <- c(inhibitor.vect, inhibitor.buffer)
        cat(inhibitor.out, file = here(inhibitors, inhibitor.out.name), sep = '\n', append = T)
        
      }
      
    }
    
  }
  
  inhibitor.log <- paste0(primer, '_inhibitors.rds')
  saveRDS(inhibitor.vect, here(logs, inhibitor.log))
  return(here(inhibitors, inhibitor.out.name))
  
}
arrester.advent <- function(df, adapter = F){
  
  region <- df[, c('gene', 'downstream_5.8s', 'its2', 'upstream_28s')]
  region$sequences <- paste0(region$downstream_5.8s, region$its2, region$upstream_28s)
  region <- region[, c(1, 5)]

  for(i in 1:nrow(region)){

    gene <- region[i, 1]
    vect <- strsplit(region[i, 2], '')[[1]]
    origin <- 1:30

    total.shifts <- length(vect) - 30
    index.shifts <- 0:total.shifts

    for(j in index.shifts){

      shifted.indices <- origin + j
      shifted.seq <- vect[shifted.indices] %>% paste(collapse = '')

      arrester.id <- paste0(gene, '_', shifted.indices[[1]])
      arrester.out <- paste0('>', arrester.id, '\n',
                              shifted.seq)

      if(i == 1 && j == 0){

        arrester.vect <- arrester.id
        cat(arrester.out, file = here(arresters, 'arresters.fasta'), sep = '\n')

      } else {

        arrester.buffer <- arrester.id
        arrester.vect <- c(arrester.vect, arrester.buffer)
        cat(arrester.out, file = here(arresters, 'arresters.fasta'), sep = '\n', append = T)

      }

    }

  }
  
  saveRDS(arrester.vect, here(logs, 'arresters.rds'))

}
compile.counts <- function(filename, overlap = T){
  
  hits <- read.delim(here(logs, filename), header = F, col.names = 'id')$id %>% 
    data.frame(id = .) %>% group_by(id) %>% summarize(hits = n())
  
  if(overlap == T ){
    
    hits %>% separate(id, c('gene', 'primer', 'position', 'overlap'), '_', F)
    
  } else {
    
    hits %>% separate(id, c('gene', 'position'), '_', F)
    
  }
  
}
best.blast <- function(df, ids, blocker){
  
  if(blocker == 'inhibitor'){
    
    keep <- setdiff(ids, df$id)
    
    if(length(keep) != 0){
      
      keep %<>% data.frame(id = .) %>% separate(id, c('gene', 'primer', 'position', 'overlap'), '_', F)
      keep$hits <- 0
      df %<>% filter(!(primer %in% keep$primer)) %>% bind_rows(keep)
      
    }
    
    df$primer.hits <- paste(df$primer, df$hits, sep = '_')
    df$primer.overlap <- paste(df$primer, df$overlap, sep = '_')
    
    # Only genes with the lowest number of hits are retained
    min.hits <- df %>% group_by(primer) %>% summarize(min.hits = min(hits))
    min.hits <- paste(min.hits$primer, min.hits$min.hits, sep = '_')
    
    df %<>% filter(primer.hits %in% min.hits)
    
    # Only primers with the greatest amount of overlap are retained, leaving the best
    max.overlap <- df %>% group_by(primer) %>% summarize(max.overlap = max(overlap))
    max.overlap <- paste(max.overlap$primer, max.overlap$max.overlap, sep = '_')
    
    filter(df, primer.overlap %in% max.overlap)
    
  } else {
    
    keep <- setdiff(ids, df$id)
    
    if(length(keep) != 0){
      
      data.frame(id = keep) %>% separate(id, c('gene', 'position'), '_', F)
      
    } else {
      
      min.hits <- min(df$hits)
      filter(df, hits <= min.hits)
      
    }
    
  }
  
}

# Clear, create, and set output directories ####
unlink('output', recursive = T)
logs <- here('output', 'logs')
figures <- here('output', 'figures')
inhibitors <- here('output', 'sequences', 'blockers', 'inhibitors')
blastdbs <- here('data', 'blastdbs')
unlink(blastdbs, recursive = T)
unlink(here('data', '5.8s_its_28s.fasta'))

dir.create(logs, recursive = T)
dir.create(figures, recursive = T)
dir.create(inhibitors, recursive = T)
dir.create(blastdbs, recursive = T)

sequences <- here('output', 'sequences')
arresters <- here(sequences, 'blockers', 'arresters')
dir.create(arresters, recursive = T)

# Set up UNITE and RDP databases ####
# Decompress
system2('gzip', args = paste('-dk', here('data', '*.gz')))

# ITS + 5.8S
db.args_5.8s <- paste(
  '-in', here('data', 'sh_general_release_dynamic_s_10.05.2021.fasta'),
  '-dbtype nucl',
  '-title "UNITE Fungal 5.8S + ITS"',
  '-out', here(blastdbs, 'unite_5.8s_its'),
  '-logfile', here(logs, 'unite_5.8s_its.txt')
)

system2('makeblastdb', args = db.args_5.8s)

# 28S
db.args_28s <- paste(
  '-in', here('data', 'current_Fungi_unaligned.fa'),
  '-dbtype nucl',
  '-title "RDP Fungal 28S"',
  '-out', here(blastdbs, 'rdp_28s'),
  '-logfile', here(logs, 'rdp_28s.txt')
)

system2('makeblastdb', args = db.args_28s)

# 5.8S + ITS + 28S
combo.args <- paste(here('data', 'current_Fungi_unaligned.fa'),
                    here('data', 'sh_general_release_dynamic_s_10.05.2021.fasta'),
                    '>',
                    here('data', '5.8s_its_28s.fasta')
                    )

system2('cat', args = combo.args)

db.args_combo <- paste(
  '-in', here('data', '5.8s_its_28s.fasta'),
  '-dbtype nucl',
  '-title "UNITE + RDP Fungal ITS + 5.8S + 28S"',
  '-out', here(blastdbs, '5.8s_its_28s'),
  '-logfile', here(logs, '5.8s_its_28s.txt')
)

system2('makeblastdb', args = db.args_combo)

# Remove non-compressed files
system2('rm', args = paste('-f', here('data', '*.fasta')))
system2('rm', args = paste('-f', here('data', '*.fa')))

# Define base blastn parameters ####
blastn.args <- paste('-perc_identity 50',
                     '-word_size 12',
                     '-dust no',
                     '-evalue 1000',
                     '-max_target_seqs 10000',
                     '-num_threads 4',
                     '-outfmt "6 qseqid"')

inhibitor.args_5.8s <- paste(blastn.args,
                             '-query', here(inhibitors, '5.8s_inhibitors.fasta'),
                             '-db', here(blastdbs, 'unite_5.8s_its'),
                             '-out', here(logs, '5.8s_inhibitors_blastn.out')
                             )


inhibitor.args_its4.fun <- paste(blastn.args,
                                 '-query', here(inhibitors, 'its4.fun_inhibitors.fasta'),
                                 '-db', here(blastdbs, 'rdp_28s'),
                                 '-out', here(logs, 'its4.fun_inhibitors_blastn.out')
                                 )

arrester.args <- paste(blastn.args,
                       '-query', here(arresters, 'arresters.fasta'),
                       '-db', here(blastdbs, '5.8s_its_28s'),
                       '-out', here(logs, 'arresters_blastn.out')
)

# Randomly generate ITS2 sequences with 5.8S and 28S priming regions based on Saccharomyces cerevisiae S288C ####
mimic <- gene.generator(adapter = 'reverse', frameshift = 3)

# Generate annealing inhibition blocking primers that overlap with priming regions ####
cat.args_5.8s <- paste(initiate.inhibitors(mimic, 'fits7'), initiate.inhibitors(mimic, 'its3.kyo1'), initiate.inhibitors(mimic, 'its3'),
                           '>',
                           here(inhibitors, '5.8s_inhibitors.fasta'))

system2('cat', args = cat.args_5.8s)

initiate.inhibitors(mimic, 'its4.fun')

# Generate elongation arrest blocking primers in the randomly generated portions of sequences ####
arrester.advent(mimic)

# Compare blockers to UNITE and RDP databases ####
# 5.8S inhibitors
system2('blastn', args = inhibitor.args_5.8s)

# 28S inhibitors
system2('blastn', args = inhibitor.args_its4.fun)

# Arresters
system2('blastn', args = arrester.args)

# Identify candidate genes and blockers ####
# Inhibitors
inhibitor.hits <- lapply(c('its4.fun_inhibitors_blastn.out', '5.8s_inhibitors_blastn.out'), compile.counts) %>% bind_rows()
inhibitor.ids <- lapply(c(here(logs, 'fits7_inhibitors.rds'),
                          here(logs, 'its3.kyo1_inhibitors.rds'),
                          here(logs, 'its3_inhibitors.rds'),
                          here(logs, 'its4.fun_inhibitors.rds')
                          ), readRDS) %>% unlist(use.names = F)

best.inhibitors <- best.blast(inhibitor.hits, inhibitor.ids, blocker = 'inhibitor')

# Arresters
arrester.hits <- compile.counts('arresters_blastn.out', overlap = F)
arrester.ids <- readRDS(here(logs, 'arresters.rds'))

best.arresters <- best.blast(arrester.hits, arrester.ids, blocker = 'arrester')

# Everything below here is junk I'm scared to throw away ####

# Define ribosomal subunit DNA sequences, derived from Saccharomyces cerevisiae S288C
# Both 5.8S and LSU are trimmed to only contain forward and reverse primers of interest

# five.8s <- 'AACTTTCAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGCGAAATGCGATACGTAATGTGAATTGCAGAATTCCGTGAATCATCGAATCTTTG'
# mid_5.8s <- 'GAAATGCGATACGTAATGTGAATTGCAGAATTCC' # 34 bp, GC content = 0.382
# downstream_5.8s <- 'AACGCACATTGCGCCCCTTGGTATTCCAGGGGGCATGCCTGTTTGAGCGTCATTT' # 55 bp, GC content = 0.545
# lsu <- 'GTTTGACCTCAAATCAGGTAGGAGTACCCGCTGAACTTAAGCATATCAATAAGCGGAGGA'
# upstream_28s <- 'GTTTGACCTCAAATCAGGTAGGAGTACCCGCTGA' # 34 bp, GC content = 0.5

# Specify the sequences of primer binding sites ####
# five.8s.fun <- 'AACTTTCAACAACGGATCTCT'
# its3.all <- 'GCATCGATGAAGAACGCAGC'
# its3 <- 'GCATCGATGAAGAACGCAGC'
# its3.kyo1 <- 'ATCGATGAAGAACGCAGC'
# its4 <- 'GCATATCAATAAGCGGAGGA'
# its4.fun <- 'ACTTAAGCATATCAATAAGCGGAGG'
# its4ngs <- 'GCATATCAATAAGCGGAGGA'
