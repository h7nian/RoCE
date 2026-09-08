# Paired working-basis comparison on fixed FACE observations, not new samples.
.source_working_basis <- function(data, config = "C1") {
  if (!is.character(config) || length(config) != 1L || is.na(config) || !config %in% c("C1","C2","C3"))
    stop("config must be C1, C2 or C3")
  if (config == "C1") return(data)
  lapply(data, function(site) {
    X <- site$X
    stopifnot(is.matrix(X), all(is.finite(X)), nrow(X) == site$n,
      identical(unname(site$Z_site), unname(cbind(X,X^2))),
      ncol(site$W_outcome) == 2L*ncol(X),
      identical(unname(site$W_outcome[,ncol(X)+seq_len(ncol(X)),drop=FALSE]),unname(X^2)))
    if (config == "C2") site$W_outcome <- X else site$Z_site <- X
    site
  })
}
