g_mean = function(x, na.rm = TRUE){
  exp(sum(log(x), na.rm = na.rm) / length(x))
}

#' Index of genes to subset on for computing local distances.
#' @param counts The n by p matrix of counts.
#' @return A vector containing the column indicies to use.
get_gene_selection <- function(counts) {
  avg <- apply(counts, 2, mean)
  index_genes <- which(avg <= quantile(avg, 0.8) & avg >= quantile(avg, 0.7))
  if(length(index_genes) == 0) {
    index_genes <- 1:ncol(counts)
  }
  return(index_genes)
}

#' Normalize using centered log ratio.
#' @param counts The n by p matrix of counts.
#' @return The normalized count matrix.
normalize_clr <- function(counts) {
  counts <- counts + 0.5 # Avoid zero expression.
  lib_sizes <- apply(counts, 1, sum) # Library size for each cell.
  counts <- counts / lib_sizes # Relative expression.
  counts <- log(counts / apply(counts, 1, g_mean)) # Centered log ratio.
  return(counts)
}

#' Get local distances among cells.
#' @param counts The n by p matrix of counts.
#' @param method The method for computing distances. See ?dist for options.
#' @return An n by n distance matrix.
get_distances <- function(counts, method = "euclidean") {
  d <- dist(counts, method = method)
  d <- as.matrix(d)
  return(d)
}
#' Kernel function applied to local distances.
#'
#' An adaptive kernel is used, which applied different bandwidths to each column.
#' As a result, the kernel is not symmetric. The weight matrix is symmetrized
#' using the pairwise minimum values.
#' @param d Either a matrix of distances.
#' @param k1 Tuning parameter.
#' @return Returns a weight matrix of transformed distances.
kernel <- function(D, k1 = NULL) {
  W <- apply(D, 2, function(x) {
    sigma <- sort(x)[k1 + 1]
    W <- exp(-x^2/sigma^2) # Use adaptive kernel.
    return(W)
  })

  # Adaptive kernel is not symmetric, so the result needs to be symmetrized.
  for(i in 2:ncol(W)) {
    for(j in 1:(i - 1)) {
      # Symmetrize using min.
      val <- min(W[i, j], W[j, i])
      W[i, j] <- val
      W[j, i] <- val
    }
  }
  return(W)
}


## Step 5: Create Markov transition matrix
get_markov <- function(s, k2 = 9) {
  # Induce sparsity in similarity matrix, if k2 provided.
  if(!is.null(k2) && !is.na(k2)) {
    s <- apply(s, 1, function(x) {
      threshold <- sort(x, decreasing = TRUE)[k2 + 1]
      x[x < threshold] <- 0
      return(x)
    })
  }
  m <- apply(s, 1, function(x) x / sum(x))
  return(m)
}

## Step 6: Perform imputation using diffusion
do_imputation <- function(counts, m) {
  m %*% counts
}


#' Obtain the Diffusion operator
#'
#' @param counts An n by p matrix of counts with cells in rows and genes in columns.
#' @param index (Optional:) an index vector for the rows of counts used to compute
#' local distance.
#' @param k1 Tuning parameter.
#' @param k2 Tuning parameter.
#' @return A diffusion object.
#' @export
#' @examples
#' df <- scDynamic::sce_example_data # Example single-cell dataset.
#' counts <- SingleCellExperiment::counts(df) # Obtain counts matrix.
#' counts <- t(counts) # Transpose so that rows correspond to cells.
#' get_diffusion(counts)
get_diffusion <- function(counts, index = NULL, k1 = 10, k2 = 20) {
  if(is.null(index) || is.na(index)) {
    index <- 1:nrow(counts)
  }
  # counts_clr <- get_clr(counts)
  # d <- get_distance(counts_clr, index)
  d <- as.matrix(dist(counts[, index], method = "canberra"))

  s <- get_similarity(d, k1 = k1)
  m <- get_markov(s, k2 = k2)
  counts <- do_imputation(counts, m)
  return(counts)
}
