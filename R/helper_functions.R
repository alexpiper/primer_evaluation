# Copyright (C) 2020 Alexander M Piper
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: alexander.piper@agriculture.vic.gov.au

# Get_primer_metrics ------------------------------------------------------

get_primer_statistics <- function(x, metrics = "all", disambiguate = TRUE){
  #Get metrics
  availmetrics <- c("base_freq","clamps",
                    "tm", "homopolymer",
                    "degeneracy", "primer_length")
  if (length (metrics)== 1 && metrics =="all"){
    metrics <- availmetrics
  }
  if(any(!metrics %in% availmetrics )){
    stop("Error, invalid metric: ", metrics[!metrics %in% availmetrics])
  }
  
  # Replace inosines
  x <- x %>% stringr::str_replace_all("I", "N")
  if(disambiguate){
    query <- unlist(DECIPHER::Disambiguate(DNAStringSet(x)))
  } else {
    query <- DNAString(x)
  }
  if("base_freq" %in% metrics){
    base_freq <- letterFrequency(query,letters="ACGT", OR=0) %>%
      prop.table() %>%
      colSums() %>%
      as.matrix() %>%
      t() %>%
      as.data.frame() %>%
      dplyr::mutate(`GC%` = (G+C))
  }
  if ("clamps" %in% metrics){ # Should ambiguous bases count?
    clamps <- data.frame(
      GC_last2 = x %>% str_extract("..$") %>% str_count("G|C"),
      GC_last5 = x %>% str_extract(".....$") %>% str_count("G|C")
    )
  }
  if ("tm" %in% metrics){
    tm <- data.frame(tm= TmCalculator::Tm_NN(x, ambiguous = disambiguate))
  }  
  if("homopolymer" %in% metrics){
    homopolymer <- data.frame(
      `poly_A` = query %>% as.character() %>% purrr::map_chr(longestConsecutive, "A") %>% max(),
      `poly_T` = query %>% as.character() %>%purrr::map_chr(longestConsecutive, "T") %>% max(),
      `poly_G` = query %>% as.character() %>% purrr::map_chr(longestConsecutive, "G") %>% max(),
      `poly_C` = query %>% as.character() %>% purrr::map_chr(longestConsecutive, "C") %>% max(),
      stringsAsFactors = FALSE
    )
  }
  if("degeneracy" %in% metrics){
    degeneracy <- data.frame(degeneracy = query %>% DECIPHER::Disambiguate() %>% length())
  } 
  if("primer_length" %in% metrics){
    primer_length <- data.frame(length = x %>% nchar())
  } 
  
  out <- metrics %>% 
    purrr::map(get, envir=sys.frame(sys.parent(0))) %>%
    bind_cols() %>%
    mutate(seq = x)
  return(out)
}

# Parralel slidenucdiag ---------------------------------------------------

slideNucDiag_para <- function (DNAbin, sppVector, width, interval = 1, cores=1) {
  
  # Define parralel nucdiag function
  nucDiag_para <- function (DNAbin, sppVector, cores=1) {
    DNAbin <- as.matrix(DNAbin)
    inform <- seg.sites(DNAbin)
    sppSeqs <- lapply(unique(sppVector), function(x) which(sppVector == x))
    # Define sitecheck
    siteCheck <- function(DNAbin, spp, inform) {
      res <- vector("logical", length=length(inform))
      for (j in 1:length(inform)) {
        site <- inform[j]
        res[j] <- as.character(DNAbin[spp, site]) %in% as.character(DNAbin[-spp, site])
        res[j] <- as.logical(sum(as.numeric(res[j])))
      }
      return(res)
    }
    li <- vector("list", length=length(sppSeqs))
    if (cores > 1){
      cl <- parallel::makeCluster(cores)
      registerDoParallel(cl)
      li <- foreach(i=1:length(sppSeqs)) %dopar% {
        siteCheck(DNAbin, sppSeqs[[i]], inform=inform)
      } 
      parallel::stopCluster(cl)
      
    } else {
      for (i in 1:length(sppSeqs)) {
        li[[i]] <- NA
        for (j in 1:length(inform)) {
          li[[i]][j] <- siteCheck(sppSeqs[[i]], inform[j])
        }
      }
    }
    out <- lapply(li, function(x) inform[which(!x)])
    names(out) <- unique(sppVector)
    return(out)
  }
  
  nd <- nucDiag_para(DNAbin, sppVector, cores=cores)
  if (interval == "codons") 
    interval <- 3
  win <- seq(1, dim(DNAbin)[2] - width, by = interval) #Get all possible windows
  mat <- matrix(NA, nrow = length(nd), ncol = length(win))
  
  #setup parallel backend to use many processors
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    registerDoParallel(cl)
    mat <- foreach(i=1:length(win), .combine = cbind) %dopar% {
      sapply(nd, function(x) length(which(x %in% win[i]:(win[i] + width))))
    }
    parallel::stopCluster(cl)
  } else{
    for (i in 1:length(win)) {
      mat[, i] <- sapply(nd, function(x) length(which(x %in% win[i]:(win[i] + width))))
    }
  }
  dimnames(mat)[[1]] <- unique(sppVector)
  dimnames(mat)[[2]] <- NULL
  return(mat)
}


# PCA BIplot --------------------------------------------------------------
PCA_biplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  plot <- ggplot(data, aes_string(x=x, y=y)) + 
    geom_point(alpha=0.5, size=3) +
    geom_text(alpha=.4, size=3, aes(label=obsnames))
  plot <- plot +  
    geom_hline(yintercept = 0, linetype=2, alpha=0.5) +  
    geom_vline(xintercept = 0, linetype=2, alpha=0.5) 
  #
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  plot
}
