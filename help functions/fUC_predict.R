trend_predict <- function(x, d, h){
    n       <- length(x)
    ar_coef <- -frac_diff(c(1, rep(0, n+h-1)), d)[-1]
    x.pred <- c(x, rep(NA, h))
    for(ii in 1:h){
        x.pred[n+ii] <- x.pred[1:(n+ii-1)]%*%(ar_coef[1:(n+ii-1)][(n+ii-1):1])
    }
    return(x.pred)
}

cycle_predict <- function(x, ar, h){
    n       <- length(x)
    ar_coef <- c(-ar, rep(0, n+h-length(ar)-1))
    x.pred <- c(x, rep(NA, h))
    for(ii in 1:h){
        x.pred[n+ii] <- x.pred[1:(n+ii-1)]%*%(ar_coef[1:(n+ii-1)][(n+ii-1):1])
    }
    return(x.pred)
}
