#'  screening longitudinal OTUs
#'
#'  @export
#'  @param ana an phyloseq object with counts data
#'  @param error_rate error rate percentage for relative abundance
#'  @param pct_threshold occurance percentage threshold
OTUscreen <- function(ana, error_rate = 0.1, pct_threshold){
  cat("range of overall counts is", range( sample_sums(ana) ),"\n")
  cat("The cutoff error rate is", error_rate, "%\n")
  cat("The cutoff occurance percentage is", pct_threshold, "%\n")

  otu <- otu_table(ana)
  otu.pct <- apply(otu_table(ana), 2, function(x) x/sum(x) )
  otu.select <- apply(otu.pct, 1, function(x) sum(x > error_rate/100) / length(x)) > pct_threshold / 100
  ana1 <- prune_taxa( names(otu.select)[otu.select == T], ana)
  cat("Number of OTUs left", nrow(otu_table(ana1)), "/", nrow(otu_table(ana)), "\n"  )
  ana1
}

#' Longitudinal correlation by samples
#'
#' @export
#' @param ana an phyloseq object with counts data
#' @param method the correlation method ("pearson","kendall","spearman")
#' @param split varialbe
#' @param nromfun the function to normalize raw counts
#' @param fill.na what value to fill if correlation is NA
tscor <- function(ana, method='kendall', splitcat='Mouse', normfun, fill.na = 0) {
  # x -> phyloseq object
  options(warn=-1)
  if (!missing(normfun)) {
    ana <- transformSampleCounts(ana, normfun)
  }
  map  <- sample_data(ana)
  otus <- otu_table(ana)
  d    <- nrow(otus)
  split.otu_table <- split.data.frame

  if (!ana@otu_table@taxa_are_rows) otus <- t(otus)

  # split OTU tables by splitcat
  otu_list <- split(t(otus), map[,splitcat])

  # pairwise correlations for each OTU sub-table
  cor_list <- lapply(otu_list, cor, method=method)
  names.cor <- names(cor_list)
  names.otu <- rownames(otus)

  ## collect correlations into multidimensional array
  cor_unlist <- unlist(cor_list)
  cor_unlist[is.na(cor_unlist)] <- fill.na
  options(warn=-1)
  cor_arr <- array(cor_unlist, dim = c(d, d, length(cor_list)),
                   dimnames = list(names.otu, names.otu, names.cor) )
}

#' Microbial Temporal Dependency Comparison
#'
#' @export
#' @param ana an phyloseq object with counts data
#' @param id.var subject (id) variable name
#' @param cov.var covariates variable name
#' @param time.var time variable name
#' @param error_rate error rate percentage for relative abundance
#' @param pct_threshold occurance percentage threshold
#' @param threshold the criteria for screening OTUs (percent of observations
#'    large than 0.1% of minimal overall counts)
#' @param method the correlation method ("pearson","kendall","spearman")
#' @param dist.type
#' @param heatmap whether to draw heatmap
#' @param classify whether to draw classifier tree
#' @example
#' #Not Run
#' load("../microdata/laura.Rdata")
#' ana  <- subset_samples(genus.count , Experiment == "Trans1" & Time >= 2 & Time <=34  )
#' map <- sample_data(ana)
#' MicroTDC(ana, id.var = "Mouse", cov.var = "Group", time.var = "Time")
MicroTDC <- function(ana, id.var, cov.var, time.var, error_rate = 0.1, pct_threshold = 20, method = "spearman", dist.type = "F", heatmap = T, classify = F){

  # Screening
  ana1 <- OTUscreen(ana, error_rate, pct_threshold)
  otu <- otu_table(ana1)
  map  <- data.frame( sample_data(ana1) )
  map.unique <- unique( map[, c(id.var, cov.var)] )
  rownames(map.unique) <- map.unique[, id.var]

  n.sample <- length(unique(map[, id.var]))

  otu.cor <- tscor(ana1, method = method, splitcat = id.var)
  dist <- outer(1:n.sample, 1:n.sample, function(x,y){
    foo <- function(x,y) norm(otu.cor[,,x] - otu.cor[,,y], type = dist.type)
    mapply(foo, x,y)
  })
  rownames(dist) <- dimnames(otu.cor)[[3]]
  colnames(dist) <- dimnames(otu.cor)[[3]]


  grp <- map.unique[rownames(dist) , cov.var ]
  # grp <- rbinom(15, size = 1, prob = 0.5)
  test <- adonis(dist ~ grp)
  pvalue <- paste( round( test$aov.tab[1,c(4,6)],3),collapse = ";")
  # print(pvalue)
  n.taxa <- nrow(otu)

  if(heatmap){
    diag(dist) <- median(dist)
    colbar = colorRampPalette(c("yellow", "red"), space = "rgb")
    heatmap.2( sqrt(dist), xlab = "Sample", ylab = "Sample", Rowv = classify, Colv = classify, dendrogram = "row",
               col = colbar, density.info = "none", trace = "none", distfun = as.dist,
               main = paste0("Cor method:", method, "Dist type:", dist.type,
                             " \nThreshold:", pct_threshold, " #OTU:", n.taxa, " pvalue:", pvalue)
    )
  }
  test
  #   dist
}


