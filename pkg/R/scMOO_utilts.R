randomized.svd <- function(A,K, q, method = 'rsvd', mkl.seed = -1) {
  out <- setNames(vector("list", 3), c("u", "d", "v"))
  if (method == 'rsvd') {
    library(rsvd)
    out <- rsvd(A,K,q=q)
  }else if (method == 'rsvd-mkl') {
    library(fastRPCA)
    fastPCAOut <- fastPCA(A, k=K, its=q, l=(K+10), seed=mkl.seed)
    out$u <- fastPCAOut$U
    out$v <- fastPCAOut$V
    out$d <- diag(fastPCAOut$S)
  }else{
    stop('Method not recognized')
  }
  return(out)
}

normalize_data <- function (A) {
  #  Simple convenience function to library and log normalize a matrix

  totalUMIPerCell <- rowSums(A);
  if (any(totalUMIPerCell == 0)) {
    toRemove <- which(totalUMIPerCell == 0)
    A <- A[-toRemove,]
    totalUMIPerCell <- totalUMIPerCell[-toRemove]
    cat(sprintf("Removed %d cells which did not express any genes\n", length(toRemove)))
  }

  A_norm <- sweep(A, 1, totalUMIPerCell, '/');
  A_norm <- A_norm * 10E3
  A_norm <- log(A_norm +1);
}

choose_k <- function (A_norm,K=100, thresh=6, noise_start=80,q=2,use.mkl=F, mkl.seed =-1) {
  #  Heuristic for choosing rank k for the low rank approximation based on
  #  statistics of the spacings between consecutive singular values. Finds
  #  the smallest singular value \sigma_i such that $\sigma_i - \sigma_{i-1}
  #  is significantly different than spacings in the tail of the singular values.
  #
  #
  # Args:
  #   A_norm: The log-transformed expression matrix of cells (rows) vs. genes (columns)
  #   K: Number of singular values to compute. Must be less than the
  #   smallest dimension of the matrix.
  #   thresh: Number of standard deviations away from the ``noise'' singular
  #   values which you consider to be signal
  #   noise_start : Index for which all smaller singular values are
  #   considered noise
  #   q : Number of additional power iterations
  #   use.mkl : Use the Intel MKL based implementation of SVD. Needs to be
  #             installed from https://github.com/KlugerLab/rpca-mkl.
  #   mkl.seed : Only relevant if use.mkl=T. Set the seed for the random
  #   generator for the Intel MKL implementation of SVD. Any number <0 will
  #   use the current timestamp. If use.mkl=F, set the seed using
  #   set.seed() function as usual.
  # Returns:
  #   A list with three items
  #       1) Chosen k
  #       2) P values of each possible k
  #       3) Singular values of the matrix A_norm

  if (K > min(dim(A_norm))) {
    stop("For an m by n matrix, K must be smaller than the min(m,n).\n")
  }
  if (noise_start >K-5) {
    stop("There need to be at least 5 singular values considered noise.\n")
  }
  noise_svals <- noise_start:K
  if (!use.mkl) {
    rsvd_out <- randomized.svd(A_norm,K,q=q)
  }else {
    rsvd_out <- randomized.svd(A_norm,K,q=q, method='rsvd-mkl', mkl.seed=mkl.seed)
  }
  diffs <- rsvd_out$d[1:(length(rsvd_out$d)-1)] - rsvd_out$d[2:length(rsvd_out$d)]
  mu <- mean(diffs[noise_svals-1])
  sigma <- sd(diffs[noise_svals-1])
  num_of_sds <- (diffs-mu)/sigma
  k <- max (which(num_of_sds > thresh))
  return (list( k=k,num_of_sds = num_of_sds,d=rsvd_out$d))
}




log.normalization = function(x, percent, preprocess.only = FALSE){

  x <- as.matrix(x)

  if (preprocess.only){

    n <- dim(x)[2]

    gene.exprs.count <- rowSums(x != 0)

    x <- x[gene.exprs.count > n * percent, ]

    return(x)

  } else {

    n <- dim(x)[2]

    gene.exprs.count <- rowSums(x != 0)

    x <- x[gene.exprs.count > n * percent, ]

    sf <- colSums(x)/median(colSums(x))

    return(log(sweep(x, 2, sf, '/')+1))

  }

}



clean.data <- function(x) {

  if (!(grepl("matrix", class(x), ignore.case = TRUE))) {

    x <- Matrix::Matrix(as.matrix(x))

    message("Converting x to matrix.")

    if (!is.numeric(x)) {

      warning("Make sure x is numeric.")

    }

  }

  np <- dim(x)

  size <- as.numeric(np[1])*as.numeric(np[2])

  if(size > 2^31-1){

    inds <- split(1:np[2], ceiling(1:np[2]/1000))

    for(i in 1:length(inds)){

      x[, inds[[i]]][x[, inds[[i]]] < 0.001] <- 0

    }

  } else {

    x[x < 0.001] <- 0

  }

  if (is.null(np) | (np[2] <= 1))

    stop("x should be a matrix with 2 or more columns")

  if (min(Matrix::colSums(x)) == 0) {

    nzerocells <- sum(Matrix::colSums(x) == 0)

    x <- x[, Matrix::colSums(x) != 0]

    message("Removing ", nzerocells, " cell(s) with zero expression.")

  }

  if (is.null(rownames(x))) {

    rownames(x) <- 1:np[1]

  }

  x

}





calc.size.factor <- function(x) {

  sf <- Matrix::colSums(x)/median(Matrix::colSums(x))

  scale.sf <- 1

  list(unname(sf), scale.sf)

}




objective <- function(W, W1, W2, W3, alpha, X, Y, AY, YB, A, B, X_rank_k, lambda1, lambda2){

  objval <- ( 1/(nrow(X)*ncol(X)) * ((1-alpha)*norm(W*(X-Y), 'F')^2 + alpha*W1*norm(X-AY, 'F')^2 + alpha*W2*norm(X-YB, 'F')^2 +

                                       alpha*W3*norm(X-t(X_rank_k), 'F')^2 +

                                       lambda1*nrow(X)*ncol(X)*sum(W1*log(W1), W2*log(W2), W3*log(W3)) +

                                       alpha*lambda2*(W1*sum(abs(A)) + W2*sum(abs(B)))) )

  return(objval)
}




regularizer_define <- function(weight_matrix ,lambda1 = 1.0, lambda2 = 1e10){
  lambda1 * k_sum(k_abs(weight_matrix), axis = c(1,2)) + lambda2 * tf$linalg$trace(tf$square(weight_matrix))
}




loss_define <- function(y_true, y_pred){
  0.5 * k_sum(k_square(y_true - y_pred), axis = c(1,2))
}



calc.lambda <- function(X_count, percent){

  X <- log.normalization(X_count, percent, preprocess.only = FALSE)

  lambda <- sd(X)

}

