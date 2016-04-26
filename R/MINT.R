#'  Screening OTUs
#'
#'Keep major taxa with a predefined error rate and percentage threshold
#'
#'
#'@param ana an phyloseq object with counts/relative abundance OTU table
#'@param error_rate error rate percentage, the default is 0.1 percent
#'@param pct_threshold occurance percentage threshold, the default percentage threshold is 20 percent
#'
#'@return an phyloseq object with major taxa
#'@export
OTUscreen <- function(ana, error_rate = 0.1, pct_threshold = 20){
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

#' OTU temporal interdependence for each subjects
#'
#'
#' @param ana a phyloseq object of counts/relative abundance data.
#' @param method an option of the correlation method ("pearson","kendall","spearman"). The default method is "spearman".
#' @param subject_var a numeric vector of subject.
#' @param fill.na a number between 0 and 1 to fill the missing value. The default value is 0.
#'
#' @return a list of temporal correlation matrix for each subject.
#'
#' @export
tscor <- function(ana, method='kendall', subject_var, fill.na = 0) {
  # x -> phyloseq object
  options(warn=-1)
  map  <- sample_data(ana)
  otus <- otu_table(ana)
  d    <- nrow(otus)
  split.otu_table <- split.data.frame

  if (!ana@otu_table@taxa_are_rows) otus <- t(otus)

  # split OTU tables by subject_var
  otu_list <- split(t(otus), map[,subject_var])

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

#' Mcrobial Interdependence Non-parametric Test (MINT)
#'
#'
#' @param otu a matrix of OTU table.
#' @param id.var a vector of subjects.
#' @param cov.var a vector of covariates.
#' @param time.var a vector of time variable.
#' @param method an option of the correlation method ("pearson","kendall","spearman"). The default method is "spearman".
#' @param dist.type the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param heatmap a logical value indicating whether to draw heatmap. The default value is TRUE.
#' @param classify a logical value indicating whether to draw classifier tree. The default value is FALSE.
#' @param fill.na a number between 0 and 1 to fill the missing value. The default value is 0.
#'
#' @return This function returns typical, but limited, output for analysis of variance (general linear models).
#' \item{aov.tab}{Typical AOV table showing sources of variation, degrees of freedom, sequential sums of squares, mean squares, F statistics, partial R-squared and P values, based on N permutations.}
#' \item{coefficients}{matrix of coefficients of the linear model, with rows representing sources of variation and columns representing species; each column represents a fit of a species abundance to the linear model. These are what you get when you fit one species to your predictors. These are NOT available if you supply the distance matrix in the formula, rather than the site x species matrix}
#' \item{coef.sites}{matrix of coefficients of the linear model, with rows representing sources of variation and columns representing sites; each column represents a fit of a sites distances (from all other sites) to the linear model. These are what you get when you fit distances of one site to your predictors.}
#' \item{f.perms}{an N by m matrix of the null F statistics for each source of variation based on N permutations of the data. The permutations can be inspected with permustats and its support functions.}
#' \item{model.matrix}{The model.matrix for the right hand side of the formula.}
#' \item{terms}{The terms component of the model.}
#'
#' @examples
#' load("data/mice.Rdata")
#' otu <- mice[, - (1:3)]
#' id.var   <- mice$id
#' cov.var  <- mice$group
#' time.var <- mice$time
#' MINT(otu, id.var, cov.var, time.var)
#'
#' @export
MINT <- function(otu, id.var, cov.var, time.var, method = "spearman", dist.type = "F", heatmap = T, classify = F, fill.na = 0){

  otu_list <- split(otu, id.var)
  d    <- ncol(otu)

  # pairwise correlations for each OTU sub-table
  cor_list <- lapply(otu_list, cor, method=method)
  n.sample <- length(cor_list)

  names.cor <- names(cor_list)
  names.otu <- colnames(otu)

  ## collect correlations into multidimensional array
  cor_unlist <- unlist(cor_list)
  cor_unlist[is.na(cor_unlist)] <- fill.na


  options(warn=-1)
  otu.cor <- array(cor_unlist, dim = c(d, d, length(cor_list)),
                   dimnames = list(names.otu, names.otu, names.cor) )

  dist <- outer(1:n.sample, 1:n.sample, function(x,y){
    foo <- function(x,y) norm(otu.cor[,,x] - otu.cor[,,y], type = dist.type)
    mapply(foo, x,y)
  })
  rownames(dist) <- dimnames(otu.cor)[[3]]
  colnames(dist) <- dimnames(otu.cor)[[3]]


  grp <- unique( data.frame(id.var , cov.var) )
  # grp <- rbinom(15, size = 1, prob = 0.5)
  test <- vegan::adonis(dist ~ cov.var, data = grp)
  pvalue <- paste( round( test$aov.tab[1,c(4,6)],3),collapse = ";")
  # print(pvalue)
  n.taxa <- ncol(otu)

  if(heatmap){
    diag(dist) <- median(dist)
    colbar = colorRampPalette(c("yellow", "red"), space = "rgb")
    gplots::heatmap.2( sqrt(dist), xlab = "Sample", ylab = "Sample", Rowv = classify, Colv = classify, dendrogram = "row",
                       col = colbar, density.info = "none", trace = "none", distfun = as.dist,
                       main = paste0("Cor method:", method, "Dist type:", dist.type,
                                    " #OTU:", n.taxa, " pvalue:", pvalue)
    )
  }
  test
}


#' Mcrobial Interdependence Non-parametric Test (MINT) using phyloseq data structure
#'
#' @export
#' @param ana an phyloseq object with counts/relative abundance data
#' @param id.var a vector of subjects.
#' @param cov.var a vector of covariates.
#' @param time.var a vector of time variable.
#' @param error_rate error rate percentage, the default is 0.1 percent.
#' @param pct_threshold occurance percentage threshold, the default percentage threshold is 20 percent.
#' @param method an option of the correlation method ("pearson","kendall","spearman"). The default method is "spearman".
#' @param dist.type the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param heatmap a logical value indicating whether to draw heatmap. The default value is TRUE.
#' @param classify a logical value indicating whether to draw classifier tree. The default value is FALSE.
#' @param fill.na a number between 0 and 1 to fill the missing value. The default value is 0.
#'
#' @return This function returns typical, but limited, output for analysis of variance (general linear models).
#' \item{aov.tab}{Typical AOV table showing sources of variation, degrees of freedom, sequential sums of squares, mean squares, F statistics, partial R-squared and P values, based on N permutations.}
#' \item{coefficients}{matrix of coefficients of the linear model, with rows representing sources of variation and columns representing species; each column represents a fit of a species abundance to the linear model. These are what you get when you fit one species to your predictors. These are NOT available if you supply the distance matrix in the formula, rather than the site x species matrix}
#' \item{coef.sites}{matrix of coefficients of the linear model, with rows representing sources of variation and columns representing sites; each column represents a fit of a sites distances (from all other sites) to the linear model. These are what you get when you fit distances of one site to your predictors.}
#' \item{f.perms}{an N by m matrix of the null F statistics for each source of variation based on N permutations of the data. The permutations can be inspected with permustats and its support functions.}
#' \item{model.matrix}{The model.matrix for the right hand side of the formula.}
#' \item{terms}{The terms component of the model.}
#'
#' @examples
#' #Not Run
#' load("../microdata/laura.Rdata")
#' ana  <- subset_samples(genus.count , Experiment == "Trans1" & Time >= 2 & Time <=34  )
#' map <- sample_data(ana)
#' MicroTDC(ana, id.var = "Mouse", cov.var = "Group", time.var = "Time")
MINT_phyloseq <- function(ana, id.var, cov.var, time.var, error.rate = 0.1, pct.threshold = 20, method = "spearman", dist.type = "F", heatmap = T, classify = F, fill.na=0){

  # Screening
  ana1 <- OTUscreen(ana, error.rate, pct.threshold)
  otu <- otu_table(ana1)
  map  <- data.frame( sample_data(ana1) )
  map.unique <- unique( map[, c(id.var, cov.var)] )
  rownames(map.unique) <- map.unique[, id.var]

  n.sample <- length(unique(map[, id.var]))

  otu.cor <- tscor(ana1, method = method, subject_var = id.var, fill.na = fill.na)
  dist <- outer(1:n.sample, 1:n.sample, function(x,y){
    foo <- function(x,y) norm(otu.cor[,,x] - otu.cor[,,y], type = dist.type)
    mapply(foo, x,y)
  })
  rownames(dist) <- dimnames(otu.cor)[[3]]
  colnames(dist) <- dimnames(otu.cor)[[3]]


  grp <- map.unique[rownames(dist) , cov.var ]
  # grp <- rbinom(15, size = 1, prob = 0.5)
  test <- vegan::adonis(dist ~ grp)
  pvalue <- paste( round( test$aov.tab[1,c(4,6)],3),collapse = ";")
  # print(pvalue)
  n.taxa <- nrow(otu)

  if(heatmap){
    diag(dist) <- median(dist)
    colbar = colorRampPalette(c("yellow", "red"), space = "rgb")
    gplots::heatmap.2( sqrt(dist), xlab = "Sample", ylab = "Sample", Rowv = classify, Colv = classify, dendrogram = "row",
               col = colbar, density.info = "none", trace = "none", distfun = as.dist,
               main = paste0("Cor method:", method, "Dist type:", dist.type,
                             " \nThreshold:", pct_threshold, " #OTU:", n.taxa, " pvalue:", pvalue)
    )
  }
  test
  #   dist
}


