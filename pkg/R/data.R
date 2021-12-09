#' Peripheral blood mononuclear cells (PBMC)

#'

#' This is the PBMC dataset sequenced using platform CEL-Seq2 (GSE132044), thus we rename it as "PBMC_CL".
#'
#' We download the preprocessed data from https://doi.org/10.5281/zenodo.3357167.
#'
#' The data contains 20,041 genes and 243 cells with true cell labels. We then select 2,000 highly variable genes using Seurat v3.2 before imputation.

#'

#' @examples

#' data("PBMC_CL")

#' @author Ke Jin, \email{kej13@mails.ccnu.edu.cn}

#' @references Ding, J. et al (2020). Systematic comparison of single-cell and single-nucleus rna-sequencing methods. \emph{Nat. Biotechnol.}, 38:737-746.

#'

"PBMC_CL"
