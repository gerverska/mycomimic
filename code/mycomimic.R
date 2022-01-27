# SYNTHETIC ITS2 GENE FRAGMENT GENERATION ####
# This gene fragment is compatible with the following fungal primers:
# fITS7, ITS3, ITS3_KYO1, 5.8S-Fun, ITS4, ITS4-Fun, ITS4ngs
# Take care to size ITS2 appropriately to your needs!

# This approach improves on J. Palmer's "amptk_synthetic_mock.py" by:
# > increasing the length of the spike-in to facilitate protocol development and downstream library re-balancing
# > generating sequences with accurate GC%
# > better sequence shuffling
# > adding blocking primer regions to allow for post-extraction tweaking of spike-in presence
# > including the option to insert adapter sequences between universal primers to allow for the resizing of synthetic amplicons at the second stage of PCR

# Load packages ####
library(tidyverse)
library(magrittr)
library(here)

# Load functions ####
# Make sure defaults are acceptable
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
seq.shuffle <- function(string, swap.rate = 0.45, shuffles = 10){
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
gene.generator <- function(perm = 200, its2.bp = 250, its2.gc = 0.5, adapter = 'none', frameshift = 'none'){
  
  if(adapter != 'none' && frameshift != 'none'){
    
    f.adapter.rc <- 'CTGTCTCTTATACACATCTGACGCTGCCGACGA' # 33 bp, GC content = 0.515
    r.adapter.rc <- 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC' # 34 bp, GC content = 0.529
    
    fs.seq <- rep('N', frameshift) %>% paste(collapse = '') 
    
    its2.part.bp <- its2.bp / 2
    its2.part.gc <- its2.gc
    
  }
  
  # Set sequence generation parameters
  # upstream_5.8s.bp <- 49 # Fixed
  # upstream_5.8s.gc <- 0.469
  mid_5.8s.bp <- 34 # Random
  mid_5.8s.gc <- 0.382
  # fits7.bp <- 19 # Fixed
  # fits7.gc <- 0.368
  downstream_5.8s.bp <- 55 # Random
  downstream_5.8s.gc <- 0.545
  # its2.bp <- 500 # Random
  # its2.gc <- 0.5
  upstream_28s.bp <- 34 # Random
  upstream_28s.gc <- 0.5
  # its4.all.bp <- 26 # Fixed
  # its4.all.gc <- 0.385
  
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
      
      its2.pt.1 <- random.seq.gc(its2.part.bp, its2.part.gc) %>%
        seq.shuffle(shuffles = 10)
      
      its2.pt.2 <- random.seq.gc(its2.part.bp, its2.part.gc) %>%
        seq.shuffle(shuffles = 10)
      
      if(adapter == 'forward'){
        
        adapter.seq <- f.adapter.rc
        
      } else {
        
        adapter.seq <- r.adapter.rc
        
      }
      
      its2 <- paste0(its2.pt.1, fs.seq, adapter.seq, its2.pt.2)
      
    } else {
      
      # adapter.seq <- 'none'
      # fs.seq <- 'none'
      # its2.part.bp <- 'none'
      # its2.part.gc <- 'none'
      its2.pt.1 <- 'none'
      its2.pt.2 <- 'none'
      
      its2 <- random.seq.gc(its2.bp, its2.gc) %>%
        seq.shuffle(shuffles = 20)
      
    }
    
    upstream_28s <- random.seq.gc(upstream_28s.bp, upstream_28s.gc) %>%
      seq.shuffle(shuffles = 10)
    
    # Combine strings to produce the full, synthetic ITS
    mimic.name <- paste0('MycoMimic.', i)
    # mimic.seq <- paste0(upstream_5.8s, mid_5.8s, fits7, downstream_5.8s, its2, upstream_28s, its4.all)
    
    # Prepare FASTA output for the full gene block, as well as separate outputs for each randomly generated component
    # mimic.out <- paste0('>', mimic.name, ';adapter=', adapter, ';frameshift=', frameshift,
    #                     ';bp_ITS2=', its2.bp, ';GC_ITS2=', its2.gc, '\n',
    #                     mimic.seq)
    # 
    # fasta.out <- paste0('mycomimic_', adapter, '.adapter.fasta')
    
    if(i == 1){
      
      mimic <- data.frame(gene = mimic.name,
                          mid_5.8s = mid_5.8s,
                          downstream_5.8s = downstream_5.8s,
                          its2.pt.1 = its2.pt.1,
                          its2.pt.2 = its2.pt.2,
                          its2 = its2,
                          upstream_28s = upstream_28s
                          )
      
      # cat(mimic.out, file = here(sequences, fasta.out), sep = '\n')
      
    } else {
      
      mimic.buffer <- data.frame(gene = mimic.name,
                                 mid_5.8s = mid_5.8s,
                                 downstream_5.8s = downstream_5.8s,
                                 its2.pt.1 = its2.pt.1,
                                 its2.pt.2 = its2.pt.2,
                                 its2 = its2,
                                 upstream_28s = upstream_28s
                                 )
      mimic <- rbind(mimic, mimic.buffer)
      
      # cat(mimic.out, file = here(sequences, fasta.out), sep = '\n', append = T)
      
    }
    
  }
  
  # Output a reference data frame with each sequence component
  return(mimic)
  
}
build.blockers <- function(df, primer){
  
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
  
  blocker.out.name <- paste0(primer, '_blockers.fasta')
  
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
      blocker.id <- paste(gene, primer, shifted.indices[[1]], current.overlap, sep = '_')
      blocker.out <- paste0('>', blocker.id, '\n',
                            shifted.seq)
      
      if(i == 1 && j == 0){
        
        blocker.vect <- blocker.id
        cat(blocker.out, file = here(sequences, blocker.out.name), sep = '\n')
        
      } else {
        
        blocker.buffer <- blocker.id
        blocker.vect <- c(blocker.vect, blocker.buffer)
        cat(blocker.out, file = here(sequences, blocker.out.name), sep = '\n', append = T)
        
      }
      
    }
    
  }
  
  blocker.log <- paste0(primer, '_blockers.rds')
  saveRDS(blocker.vect, here(logs, blocker.log))
  return(here(sequences, blocker.out.name))
  
}
compile.counts <- function(filename, overlap = T){
  
  hits <- read.delim(here(logs, filename), header = F, col.names = 'id')$id %>% 
    data.frame(id = .) %>% group_by(id) %>% summarize(hits = n())
    
  hits %<>% separate(id, c('gene', 'primer', 'position', 'overlap'), '_', F)
  
}
best.blast <- function(df, ids){
    
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
  
  min.max <- filter(df, primer.overlap %in% max.overlap) %>% group_by(primer) %>%
    summarize(id = sample(id, 1))
  min.max$gene <- str_extract(min.max$id, '^[[:alpha:]]+.[[:digit:]]+')
  
  min.max
  
}

# Clear, create, and set output directories ####
unlink('output', recursive = T)
logs <- here('output', 'logs')
seq.out <- here('output', 'sequences')
seq.in <- here('input', 'sequences')
blastdbs <- here('output', 'blastdbs')
unlink(blastdbs, recursive = T)

dir.create(logs, recursive = T)
dir.create(seq.out, recursive = T)
dir.create(blastdbs, recursive = T)

# Set up UNITE and RDP databases ####
# Decompress
system2('gzip', args = paste('-dk', here(seq.in, '*.gz')))

# ITS + 5.8S
db.args_5.8s <- paste(
  '-in', here(seq.in, 'sh_general_release_dynamic_s_10.05.2021.fasta'),
  '-dbtype nucl',
  '-title "UNITE Fungal 5.8S + ITS"',
  '-out', here(blastdbs, 'unite_5.8s_its'),
  '-logfile', here(logs, 'unite_5.8s_its.txt')
)

system2('makeblastdb', args = db.args_5.8s)

# 28S
db.args_28s <- paste(
  '-in', here(seq.in, 'current_Fungi_unaligned.fa'),
  '-dbtype nucl',
  '-title "RDP Fungal 28S"',
  '-out', here(blastdbs, 'rdp_28s'),
  '-logfile', here(logs, 'rdp_28s.txt')
)

system2('makeblastdb', args = db.args_28s)

# 5.8S + ITS + 28S
# combo.args <- paste(here('data', 'current_Fungi_unaligned.fa'),
#                     here('data', 'sh_general_release_dynamic_s_10.05.2021.fasta'),
#                     '>',
#                     here('data', '5.8s_its_28s.fasta')
#                     )
# 
# system2('cat', args = combo.args)
# 
# db.args_combo <- paste(
#   '-in', here('data', '5.8s_its_28s.fasta'),
#   '-dbtype nucl',
#   '-title "UNITE + RDP Fungal ITS + 5.8S + 28S"',
#   '-out', here(blastdbs, '5.8s_its_28s'),
#   '-logfile', here(logs, '5.8s_its_28s.txt')
# )
# 
# system2('makeblastdb', args = db.args_combo)

# Remove non-compressed files
unlink(here(seq.in, '*.fasta'))
unlink(here(seq.in, '*.fa'))

# Define blastn parameters ####
blastn.args <- paste('-perc_identity 50',
                     '-word_size 12',
                     '-dust no',
                     '-evalue 1000',
                     '-max_target_seqs 10000',
                     '-num_threads 4',
                     '-outfmt "6 qseqid"')

# For ITS3, ITS3-KYO1, and fITS7 annealing inhibition blocking primers
blocker.args_5.8s <- paste(blastn.args,
                             '-query', here(seq.out, '5.8s_blockers.fasta'),
                             '-db', here(blastdbs, 'unite_5.8s_its'),
                             '-out', here(logs, '5.8s_blockers_blastn.out')
                             )

# For ITS4-Fun annealing inhibition blocking primers
blocker.args_its4.fun <- paste(blastn.args,
                                 '-query', here(seq.out, 'its4.fun_blockers.fasta'),
                                 '-db', here(blastdbs, 'rdp_28s'),
                                 '-out', here(logs, 'its4.fun_blockers_blastn.out')
                                 )

# Randomly generate ITS2 sequences with 5.8S and 28S priming regions based on Saccharomyces cerevisiae S288C ####
# Adjust adapter and frameshift arguments as needed
mimic <- gene.generator()

# Define ribosomal subunit DNA sequences, derived from Saccharomyces cerevisiae S288C
# Both 5.8S and LSU are respectively trimmed upstream and downstream to only contain forward and reverse primers of interest

upstream_5.8s <- 'AACTTTCAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGC' # 49 bp, GC content = 0.469
fits7 <- 'GTGAATCATCGAATCTTTG' # 19 bp, GC content = 0.368
# fits7 <- 'GTGARTCATCGAATCTTTG'
its4.all <- 'ACTTAAGCATATCAATAAGCGGAGGA' # 26 bp, GC content = 0.385

mimic$upstream_5.8s <- upstream_5.8s
mimic$fits7 <- fits7
mimic$its4.all <- its4.all

# Generate annealing inhibition blocking primers that overlap with priming regions ####
cat.args_5.8s <- paste(build.blockers(mimic, 'fits7'), build.blockers(mimic, 'its3.kyo1'), build.blockers(mimic, 'its3'),
                           '>',
                           here(seq.out, '5.8s_blockers.fasta'))

system2('cat', args = cat.args_5.8s)

build.blockers(mimic, 'its4.fun')

# Restriction enzyme survey ####
re <- readRDS(here('data', 're.sequences.rds'))

fastafy <- function(df, filename){
  
  for(i in 1:nrow(df)){
    
    line.out <- paste0('>', df[, 're'], '\n',
                       df[, 'sequence'])
    
    file.out <- paste0(filename, '.fasta')
    
    if(i == 1){
      
      cat(line.out, file = here(data, file.out), sep = '\n')
      
    } else {
      
      cat(line.out, file = here(seq.out, file.out), sep = '\n', append = T)
      
    }
    
  }
  
}

fastafy(re, 're.sequences.fasta')

# Compare blockers to UNITE and RDP databases ####
# 5.8S blockers
system2('blastn', args = blocker.args_5.8s)

# 28S blockers
system2('blastn', args = blocker.args_its4.fun)

# Identify candidate genes and blockers ####
blocker.hits <- lapply(c('its4.fun_blockers_blastn.out', '5.8s_blockers_blastn.out'), compile.counts) %>% bind_rows()
blocker.ids <- lapply(c(here(logs, 'fits7_blockers.rds'),
                          here(logs, 'its3.kyo1_blockers.rds'),
                          here(logs, 'its3_blockers.rds'),
                          here(logs, 'its4.fun_blockers.rds')
                          ), readRDS) %>% unlist(use.names = F)

best <- best.blast(blocker.hits, blocker.ids) %>% left_join(mimic, by = 'gene')
saveRDS(best, here(logs, 'best.mycomimic.rds'))

# Combine the best into one ####
m13f <- 'GTAAAACGACGGCCAGT'
m13r <- strsplit('CAGGAAACAGCTATGAC', '')[[1]] %>% rev() %>% paste(collapse = '') %>%
  chartr('ATGC','TACG', .)

paste0('>MycoMimic', '\n',
       m13f,
       upstream_5.8s,
       best[best$primer == 'its3', ]$mid_5.8s,
       fits7,
       best[best$primer == 'fits7', ]$downstream_5.8s,
       best[best$primer == sample(best$primer, 1), ]$its2,
       best[best$primer == 'its4.fun', ]$upstream_28s,
       its4.all,
       m13r
       ) %>% cat(., file = here(seq.out, 'mycomimic.fasta'))

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
# five.8s.fun <- 'AACTTTYRRCAAYGGATCWCT'
# its3.all <- 'GCATCGATGAAGAACGCAGC'
# its3 <- 'GCATCGATGAAGAACGCAGC'
# its3.kyo1 <- 'ATCGATGAAGAACGCAG'
# its3.kyo1 <- 'AHCGATGAAGAACRYAG'
# its4 <- 'GCATATCAATAAGCGGAGGA'
# its4.fun <- 'ACTTAAGCATATCAATAAGCGGAGG'
# its4.fun <- 'AYTTAAGCATATCAATAAGCGGAGG'
# its4ngs <- 'GCATATCAATAAGCGGAGGA'
# its4ngs <- 'GCATATCAATAAGCGSAGGA'
