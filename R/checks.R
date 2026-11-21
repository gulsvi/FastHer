###################################################################
###              Consistency of SNP names and order             ###
###################################################################
 
# Compares SNP names when they are available: names of z, row and column
# names of ld, names of the weights. Inputs without names are not compared.
.check_snp_names <- function(z, ld, weights = NULL) {
    zn <- names(z)
    ldn <- colnames(ld)
    ldr <- rownames(ld)
    if (!is.null(ldn) && !is.null(ldr) && !identical(ldn, ldr)) {
        stop("Row names and column names of 'ld' differ.", call. = FALSE)
    }
    if (!is.null(zn) && anyDuplicated(zn) > 0L) {
        stop("'z' has duplicated SNP names.", call. = FALSE)
    }
    if (!is.null(zn) && !is.null(ldn) && !identical(zn, ldn)) {
        stop("SNP names of 'z' and 'ld' do not match: ", sum(!zn %in% ldn),
             " SNPs of 'z' are not in 'ld'. If the SNPs are the same, their order differs; ",
             "put 'z' and the rows and columns of 'ld' in the same order.", call. = FALSE)
    }
    wn <- names(weights)
    if (!is.null(zn) && !is.null(wn) && !identical(zn, wn)) {
        stop("SNP names of 'weights' and 'z' do not match.", call. = FALSE)
    }
    invisible(TRUE)
}
