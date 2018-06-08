#' diss.relative
#'
#' @param D 
#'
#' @return estimates a relative distance matrix based on appendix S1 (Ypma et al. 2013)
#' @export
diss.relative <- function(D){
        d <- matrix(0, nrow = nrow(D), ncol = ncol(D))
        for(i in 1:(nrow(D)-1)){
                for(j in (i+1):nrow(D)){
                        if(D[i, j] != 0) {
                                x <- sum(D[i,] < D[i,j] & D[,j] < D[j,i])
                                x2 <- (sum(D[i,] == 0) + sum(D[,j] == 0))/2
                                ties.i <- sum(D[i,] < D[i,j] & D[i,] !=0 & D[,j] == D[j,i])/2
                                ties.j <- sum(D[,j] < D[i,j] & D[,j] !=0 & D[i,] == D[j,i])/2
                                ties <- sum(D[i,] == D[i,j] & D[,j] == D[j,i])/3
                                d[i,j] <- x + x2 + ties.i + ties.j + ties
                        } else {
                                d[i,j] <- (sum(D[i,] == 0) + 1)/3
                        }
                        d[j,i] <- d[i,j]
                }
        }
        return(d)
}