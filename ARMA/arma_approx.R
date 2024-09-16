library(CFFpack)

setwd("/Users/tobias/Dokumente/Projekte/Filtering unknown persistence/R/code/ARMA")
CFFpack::ARMA_approx(p=3, q=3, n=100, ncl=4)
CFFpack::ARMA_approx(p=3, q=3, n=200, ncl=4)
CFFpack::ARMA_approx(p=3, q=3, n=300, ncl=4)
