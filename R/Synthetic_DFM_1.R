#====Synthetic DFM with idio structures====
#fixed
n = 100 # number of observations
N = 60 # number of variables 
p = 1 # VAR(p) model for factors 
s = 1 # weight of true non-zero loadings

#tweaking parameters
r = 10 # number of factors 
rho = 0.7 # VAR(p) auto-correlation coefficient  (0.3,0.8)
d = 0.8 # serial correlation for idio (0.3,0.8)
tau = 0.8 # cross correlation for idio (0.3,0.8)

x <- 1
repeat {
  A = rho*diag(r)
  Sig_u = diag(r) #fixed sig_u = 1
  factors = tsDyn::VAR.sim( B = A, n = n, include = "none",  varcov = Sig_u)
  
  Lambda_dense = matrix(rnorm(N*r), ncol = r)
   
  D = d*diag(N)
  Ip <- diag(N)
  Tau = tau*Ip
  
  e = tsDyn::VAR.sim( B = D, n = n, include = "none",  varcov = Tau)
  
  X = factors %*% t(Lambda_dense) + e
  
  form = sprintf('subject_%s.csv', x)
  write.csv(X, file = form, row.names = FALSE)
  x = x+1
  if (x == 101) {
    break
  }
}
