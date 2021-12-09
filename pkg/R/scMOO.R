#' use scMOO to impute dropout values in scRNA-seq data
#'
#' @param Y_count An expression count matrix. The rows correspond to genes and
#' the columns correspond to cells. Can be sparse.
#'
#' @param W W_{gc} is set to the ratio of the number of non-zero elements in Y to the number
#' of zero elements in Y if Y_{gc} = 0, and to 1 otherwise.
#'
#' @param k The rank of the low-rank approximation matrix.
#'
#' @param q The parameter for number of additional power iterations.
#'
#' @param percent The expression count matrix is preprocessed by filtering out the genes
#' expressed in at most percent*\eqn{100\%} of the cells.
#'
#' @param lambda1 Tuning parameter for entropy regularizer.
#'
#' @param lambda2 Tuning parameter to facilitate feature selection and regularization.
#'
#' @param lambda3 Tuning parameter to penalize the diagonal elements of the parameter to
#' eliminate the trivial solution of representing an expression level as a linear combination
#' of itself.
#'
#' @param alpha Tuning parameter to balance the error between the imputed and
#' observed data and the error of model fitting data.
#'
#' @param MAX_ITER Maximum iteration of scMOO.
#'
#' @param ABSTOL Absolute tolerance of circulation.
#'
#' @param learning_rate A hyper-parameter that controls the speed of adjusting the weights of the network
#' with respect to the loss gradient.
#'
#' @param epochs The number of the entire training set going through the entire network.
#'
#' @param verbose Whether to output the value of metrics at the end of each epoch. Default is TRUE.
#'
#'
#'
#' @return If `estimates.only = TRUE', then a matrix of scMOO estimates.
#'
#' If `estimates.only = FALSE', a list with the following components
#'
#' \item{\code{estimate}}{Recovered (normalized) expression.}
#'
#' \item{\code{size.factor}}{Size factor used for normalization.}
#'
#' \item{\code{pred.time}}{Total time for scMOO estimation.}
#'
#' \item{\code{alpha}}{Tuning parameter to balance the error between the imputed and
#' observed data and the error of model fitting data.}
#'
#' \item{\code{w}}{The combination weights W1, W2 and W3.}
#'
#'
#' @export
#'
#'
#' @import keras
#'
#' @import tensorflow
#'
#' @import rsvd
#'
#' @author Ke Jin, \email{kej13@mails.ccnu.edu.cn}
#'
#' @examples
#'
#' data("PBMC_CL")
#'
#' result <- scMOO(PBMC_CL, percent = 0, estimates.only = TRUE)
#'
scMOO <- function(Y.count, W=NULL, K=0, q=10, percent = 0.05,

                lambda1=NULL, lambda2=NULL, lambda3=1e10, alpha=NULL,

                MAX.ITER = 4, ABSTOL = 1e-3, learning.rate = 0.0001,

                epochs = 100, verbose = TRUE, estimates.only = FALSE){


  # calculating the tunning parameter lambda

  if (is.null(lambda2)) {

    message("Calculating the penalty parameter in the lasso models ...")

    lambda2 <- calc.lambda(Y.count, percent)

    message("Done!")

  }


  # preprocessing and log-normalization

  message("Starting preprocessing and log-normalization ...")

  Y <- log.normalization(Y.count, percent, preprocess.only = FALSE)

  Y.count <- log.normalization(Y.count, percent, preprocess.only = TRUE)

  message("Done!")


  # compute weights
  if (is.null(W)) {

    m = nrow(Y)
    n = ncol(Y)

    W <- matrix(1, m, n)
    rate <- round(sum(Y!=0)/sum(Y==0), 2)

    W[Y==0] <- rate
  }


  Y.count <- clean.data(Y.count)

  ngenes <- nrow(Y.count)

  ncells <- ncol(Y.count)

  gene.names <- rownames(Y.count)

  cell.names <- colnames(Y.count)


  # assign size factor

  sf.out <- calc.size.factor(Y.count)

  sf <- sf.out[[1]]

  scale.sf <- sf.out[[2]]



  result <- list()

  result$size.factor <- scale.sf*sf



  message("Imputation starts ...")

  message("iter", ' objective')

  pred.st <- Sys.time()

  k <- 1

  history.objval <- c()


  # initialize X
  X <- Y


  batch.size <- nrow(Y)

  B <- keras_lasso_regression(X, Y, epochs = epochs, batch_size = batch.size, lambda1=lambda2, lambda2=lambda3,

                              learning_rate = learning.rate, verbose = verbose)

  batch.size <- ncol(Y)

  A <- t(keras_lasso_regression(t(X), t(Y), epochs = epochs, batch_size = batch.size, lambda1 = lambda2, lambda2 = lambda3,

                                learning_rate = learning.rate, verbose = verbose))

  if ( K ==0 ) {
    k_choice <- choose_k(t(X))
    K <-  k_choice$k
    cat(sprintf("Chose K=%d\n",K))
  }

  cat("Randomized SVD\n")
  fastDecomp_noc <- randomized.svd(t(X),K,q=q)

  X_rank_k <- fastDecomp_noc$u[,1:K]%*%diag(fastDecomp_noc$d[1:K])%*% t(fastDecomp_noc$v[,1:K])

  X_rank_k[X_rank_k<0]=0

  ########
  AY <- A%*%Y
  YB <- Y%*%B


  X <- ( 0.5*W^2*Y + 0.5*(1/3*pmax(AY,0)+1/3*pmax(YB,0)+1/3*t(X_rank_k)) )/ (0.5*W^2+0.5)


  message("Calculating the tuning parameter alpha ...")
  if (is.null(alpha)) {

    alpha <- norm(W*(X-Y), 'F')^2/(1/3*norm(X-AY, 'F')^2 + 1/3*norm(X-YB, 'F')^2
                                   + 1/3*norm(X-t(X_rank_k), 'F')^2 + norm(W*(X-Y), 'F')^2)
  }


  if (is.null(lambda1)) {

    message("Calculating the tuning parameter for entropy regularizer ...")

    lambda1 <- alpha * mean( norm(X-AY, 'F')^2, norm(X-YB, 'F')^2,
                             norm(X-t(X_rank_k), 'F')^2)/(ngenes*ncells)

    message("Done!")

  }


  while (k <= MAX.ITER){


    batch.size <- nrow(Y)


    # using keras with SGD algorithm to calculate B or A

    B <- keras_lasso_regression(Y=X, X=Y, epochs = epochs, batch_size = batch.size, lambda1=lambda2, lambda2=lambda3,

                                learning_rate = learning.rate, verbose = verbose)


    batch.size = ncol(Y)

    A <- t(keras_lasso_regression(Y=t(X), X=t(Y), epochs = epochs, batch_size = batch.size, lambda1 = lambda2, lambda2 = lambda3,

                                  learning_rate = learning.rate, verbose = verbose))


    # calculating the low-rank approximation
    k_choice <- choose_k(t(X))
    K <-  k_choice$k
    cat(sprintf("Chose K=%d\n",K))


    cat("Randomized SVD\n")
    fastDecomp_noc <- randomized.svd(t(X),K,q=q)

    X_rank_k <- fastDecomp_noc$u[,1:K]%*%diag(fastDecomp_noc$d[1:K])%*% t(fastDecomp_noc$v[,1:K])

    X_rank_k[X_rank_k<0]=0


    # compute W1, W2 ,W3
    AY <- A%*%Y
    YB <- Y%*%B

    N1 <- alpha*norm(X-AY, 'F')^2/(lambda1*ngenes*ncells)
    N2 <- alpha*norm(X-YB, 'F')^2/(lambda1*ngenes*ncells)
    N3 <- alpha*norm(X-t(X_rank_k), 'F')^2/(lambda1*ngenes*ncells)

    W_sum <- sum(exp(-N1), exp(-N2), exp(-N3))
    W1 <- exp(-N1)/W_sum
    W2 <- exp(-N2)/W_sum
    W3 <- exp(-N3)/W_sum
    cat(sprintf("W1=%.4f\n",W1))
    cat(sprintf("W2=%.4f\n",W2))
    cat(sprintf("W3=%.4f\n",W3))


    X <- ( (1-alpha)*W^2*Y + alpha*(W1*AY + W2*YB + W3*t(X_rank_k)) )/ ((1-alpha)*W^2 + alpha)


    history.objval[k] <- objective(W=W, W1=W1, W2=W2, W3=W3, alpha=alpha, X=X, Y=Y, AY=AY, YB=YB,
                                    A=A, B=B, X_rank_k=X_rank_k, lambda1=lambda1, lambda2=lambda2)


    message(k, '    ', round(history.objval[k], 3))

    if (k > 1 && abs(history.objval[k] - history.objval[k-1]) < ABSTOL) break

    k <- k + 1



  }



  Xhat <- pmax(X, 0)



  obj.prior <- history.objval


  pred.time <- Sys.time() - pred.st


  result$estimate <- exp(Xhat)-1


  result$pred.time <- pred.time

  result$alpha <- alpha

  result$w <- c(W1, W2, W3)

  message("Done!")

  message("Finish time: ", Sys.time())

  message("Total time: ", format(result$pred.time))

  if (!estimates.only) {

    class(result) <- "scMOO"

    result

  } else {

    result$estimate

  }

  class(result) <- "scMOO"

  result

}
