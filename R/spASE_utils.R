expit <- function(x) {
  expx <- exp(x)
  return(expx/(1+expx))
}

dexpit <- function(x) {
  expitx <- expit(x)
  return(expitx*(1 - expitx))
}

texpit <- function(x) {
  expitx <- expit(x)
  return(expitx*(1 - expitx)*(1 - 2*expitx))
}
