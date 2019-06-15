#' different way of integration
#' This function generates the values of integration.
#' @param f  a function which you want to integrate.
#' @param L  lower bound.
#' @param U  upper bound.
#' @param type  different way of integration "simpson"(default),"trapezoid","riemann"
#' @param n  cutting into n pieces.
#' @param t  simulating t times.
#' @return a scalar equivalent to integrate from lower bound to upper bound.
#' @export
#' @seealso \code{integrate}
#' @references
#' http://www.tf.uns.ac.rs/~omorr/radovan_omorjan_003_prII/r-examples/spuRs/spuRs-2ed.pdf
#' @examples
#' ## integrate from 1 to 2 by simpson method
#' f=function(x) 4*x^3
#' diffintegrate(f,1,2)
#'
#'## integrate by trapezoid method
#' f=function(x) 4*x^3
#' diffintegrate(f,1,2,type="trapezoid",n=50,t=100)


diffintegrate <- function(f, L, U, type="simpson",n=500, t=50) {

  if(type=="trapezoid"){
    u=matrix(0,t,n+2)
    tt=NULL
    k=matrix(0,t,n)
    y=NULL
    for (j in 1:t){
      u[j,] <- c(L,sort(runif(n,L,U)),U)
      for (i in 1:n+1) {
        tt=c(tt,((f(u[j,(i+1)])+f(u[j,i]))*(u[j,i+1]-u[j,i]))/2)
      }
      y=cbind(y,sum(tt[(1+(j-1)*n):(n+n*(j-1))]))
    }
    return(mean(y))
  }else if(type=="simpson"){
    sum_A <- rep(NA,t)
    A <- rep(NA,t*n)
    n1 <- n-1
    h <- rep(NA,t*n1)
    for (i in 1:t) {
      h[seq((n1*i)-(n1-1),n1*i)] <- sort(runif(n1, min = L, max = U))
      a <- c(L,h[seq((n1*i)-(n1-1),n1*i)])
      b <- c(h[seq((n1*i)-(n1-1),n1*i)],U)
      A[seq((n*i)-(n-1),n*i)] <- (b-a)/6*(f(a)+4*f((a+b)/2)+f(b))
      sum_A[i] <- sum(A[seq((n*i)-(n-1),n*i)])
    }
    r <- mean(sum_A)
    return(r)
  }else if(type=="riemann"){
    A <- rep(NA,t)
    for (j in 1:t) {
      n1 <- (n-1)
      h_upper <- rep(NA,n)
      h_lower <- rep(NA,n)
      h <- c(L, sort(runif(n1, min = L, max = U)), U)
      y <- rep(NA,n+1)
      for (q in 1:(n+1)) {
        y[q] <- f(h[q])
      }
      for (i in 1:n) {
        h_upper[i] <- max(y[i],y[i+1])
        h_lower[i] <- min(y[i], y[i+1])
        dis[i] <- h[i+1]-h[i]
      }
      dis <- round(dis, digits = 6)
      h_upper <- round(h_upper, digits = 6)
      h_lower <- round(h_lower, digits = 6)
      dis_h_u <- rep(NA,n)
      dis_h_l <- rep(NA,n)
      for (i in 1:n) {
        dis_h_u[i] <- sum(dis[i]*h_upper[i])
        dis_h_l[i] <- sum(dis[i]*h_lower[i])
      }
      A_upper <- sum(dis_h_u)
      A_lower <- sum(dis_h_l)
      A[j] <- mean(A_upper, A_lower)
    }
    r <- mean(A)
    return(r)
  }

}


