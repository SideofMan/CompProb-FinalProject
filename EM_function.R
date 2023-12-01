EM_function <- function(X, m, n = 50){
  # performs the Expectation-maximization algorithm on data X
  # with m mixing components, and n = 50 iterations at most
  
  N = length(X)
  
  # initialize the parameters
  
  pi = matrix(0,n,m); mu = matrix(0,n,m); sigma = matrix(0,n,m)
  pi[1,] = rep(1/m,m)
  mu[1,] = seq(from = min(X), to = max(X), length.out = m) # space them out evenly
  sigma[1,] = rep(diff(range(X))/(6*m),m) # spread them out evenly
  
  theta = list(pi = pi, mu = mu, sigma = sigma)
  
  iter = 1
  
  while (iter < n){
    # calculate the conditional probabilities
    p = matrix(0,N,m)
    for (j in 1:m){
      p[,j] = theta$pi[iter,j]*dnorm(X, mean = theta$mu[iter,j], sd = theta$sigma[iter,j])
    }
    
    p.hat = t(apply(p, 1, function(x){x/sum(x)}))
    
    # calculate the new parameters
    
    new_pi = apply(p.hat,2,sum)/N
    new_mu = apply(p.hat,2,function(x){sum(x*X)/sum(x)})
    new_sigma = sqrt(apply(p.hat,2,function(x){sum(x*(X-sum(x*X)/sum(x))^2)/sum(x)}))
    
    theta$pi[iter+1,]=new_pi
    theta$mu[iter+1,]=new_mu
    theta$sigma[iter+1,]=new_sigma
    
    if (iter == 1){
      LogLikelihood = c(sum(log(rowSums(p))))
    } else {
      LogLikelihood = c(LogLikelihood, sum(log(rowSums(p))))
    }
    
    if (iter > 1 & abs((LogLikelihood[iter+1] - LogLikelihood[iter])/LogLikelihood[iter]) < 1e-5){
      break
    }
    
    iter = iter + 1
  }
  
  return(theta)
}