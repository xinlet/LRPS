#####Synthetic DFM without idio structures######
#fixed
n = 100 # number of observations
N = 60 # number of variables 
p = 1 # VAR(p) model for factors 
s = 1 # weight of true non-zero loadings

#tweaking parameters
r = 10 # number of factors (5,10)
rho = 0.3 # VAR(p) auto-correlation coefficient  (0.3,0.8)

x <- 1
repeat {
  A = rho*diag(r)
  Sig_u = diag(r) #fixed sig_u = 1
  factors = tsDyn::VAR.sim( B = A, n = n, include = "none",  varcov = Sig_u)
  
  Lambda_dense = matrix(rnorm(N*r), ncol = r)

  e = matrix(rnorm(n*N), ncol = N)
 
  X = factors %*% t(Lambda_dense) + e
  
  form = sprintf('subject_%s.csv', x)
  write.csv(X, file = form, row.names = FALSE)
  x = x+1
  if (x == 101) {
    break
  }
}
