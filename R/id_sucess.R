sapply(c("rentrez", "bold", "taxize","taxizedb", "usethis", "tidyverse", "spider", "insect", "ape", "DECIPHER", "ggpubr", "RColorBrewer", "plotly", "ggforce", "seqinr", "shortread", "patchwork", "viridis","ggridges","UpSetR"), require, character.only = TRUE)


primers <- read_csv("primer_evaluation/primer_candidates.csv")
Alignments <- sort(list.files("primer_evaluation/pestfamilies/", pattern = ".fa.gz", full.names = TRUE)) # Read fasta filenames

names <- basename(Alignments) %>%
  str_replace(".fa.gz","")


# Filter to only those within alignment size & HiSeq sized
dat.passed <- primers %>%
  filter(F.Start > 0 & R.Stop < 661) %>%
  filter(amplicon < 240) %>% 
  mutate(F.seq = str_replace_all(F.seq, "I","N"))%>% #Replace Inosines with N
  mutate(R.seq = str_replace_all(R.seq, "I","N"))

p <- 1
i <- 1
dat <- vector("list", length=nrow(dat.passed))
prime <- vector("list", length=nrow(dat.passed))
dir.create("primer_evaluation/amplicons")


for (i in 1:length(Alignments)) {
  
  name <- names[i]
  seqs <- readFASTA(Alignments[i])
  for (p in 1:nrow(dat.passed)) {
    
    amplicon <- virtualPCR(seqs, up = dat.passed$F.seq[p], dat.passed$R.seq[p], rcdown = TRUE, trimprimers = TRUE, quiet = TRUE)
    if (length(amplicon) > 3) {
      
      #Filter to median - Some amplicons have primer slippage?
      seqLength <- sapply(amplicon, length)
      
      # Get most frequent value
      uniqx <- unique(na.omit(seqLength))
      if (length(uniqx) > 1 ) {message(dat.passed$F.Name[p]," and ",  dat.passed$R.Name[p], " have primer slippage for ", name)}
      freqlen <- uniqx[which.max(tabulate(match(seqLength, uniqx)))]
      
      amplicon <- as.matrix(amplicon[which(seqLength == freqlen)])
      
      # Genus and species names
      aa <- Biostrings::strsplit(dimnames(amplicon)[[1]], split = ";")
      Genus <- sapply(aa, function(x) paste(x[7], sep = "_"))
      Spp <- sapply(aa, function(x) paste(x[8], sep = "_"))
      
      Dist <- dist.dna(amplicon, pairwise.deletion = TRUE)
      
      if (any(Dist > 0)){
        closematch <- tibble(
          query = labels(amplicon),
          Spp = Spp,
          result = bestCloseMatch(Dist, Spp))  %>% 
          left_join(enframe(bestCloseMatch(Dist, Spp, names = TRUE)) %>%
                      set_colnames(c("query", "names"))%>% 
                      mutate(names = map(names, ~set_names(., paste0("closematch_",seq_along(.))))) ,
                    by="query") %>%
          unnest_wider(col=names)
        
        if (length(unique(Spp)) > 3) {
          Tr <- nj(Dist)
          maxInt <- max(Tr$edge.length[Tr$edge[, 2] > length(Tr$tip.label)])
          nodeRoot <- Tr$edge[which(Tr$edge.length == maxInt), 2]
          TrRoot <- root(Tr, node = nodeRoot, resolve.root = TRUE)
          TrRoot$tip.label <- Spp
          mono <- monophyly(TrRoot, Spp, singletonsMono = TRUE)
          
          prime[[p]] <-  tibble(
            primer = dat.passed$Name[p],
            Spp = Spp,
            query = dimnames(amplicon)[[1]],
            nn = nearNeighbour(Dist, Spp),
            nn_spp = nearNeighbour(Dist, Spp, names = TRUE),
            mono=mono[match(Spp, unique(Spp))]
          ) %>%
            left_join(closematch, by=c("query", "Spp"))
          
        } else if (length(unique(Spp)) <3) {
          prime[[p]] <- tibble(
            primer = dat.passed$Name[p],
            Spp = Spp,
            query = dimnames(amplicon)[[1]],
            nn = nearNeighbour(Dist, Spp),
            nn_spp = nearNeighbour(Dist, Spp, names = TRUE)
          ) %>%
            left_join(closematch, by=c("query", "Spp"))
        } 
      } 
    } else next()
  }
  out <- bind_rows(prime)
  write_csv(out, paste0("primer_evaluation/amplicons/",name,".csv"))
  
  dat[[i]] <- out
}