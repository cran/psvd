
#############################################################################
calcSVD <- function(X, r=2, eta=10, itmax=200, err=1e-8, mySeed=50) {
    d <- ncol(X)

    # Initial random orthonormal matrix
    set.seed(mySeed)
    w <- c(rnorm(d*r, 0,1))
    res <- mGSc(w, d, r)
    wp <- matrix(res$wp, ncol = r)

    # power 2 symmetric matrix to use in the iterative search
    t1 <- diag(c(rep(1,d))) + eta*t(X) %*% X
    tt <- t1 %*% t1

    # calculation of the eigenvectors matrix
    res <- eigenVc(tt, wp, d, r, itmax, err)
    v <- matrix(res$wc, ncol = r)
    t <- res$iter

    # results extraction for the SVD
    xw <- X %*% v
    tt <- diag(t(xw) %*% xw)
    d <- sqrt(tt)
    u <- xw %*% diag(1/d)
    list(d=d,u=u,v=v,iter=t)
}
#############################################################################
