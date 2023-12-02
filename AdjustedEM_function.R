  #I wanted to show the EM funtion in our presentation but the original version was a little lengthy, this does the exact same thing just a little more condensed 
  #So I can fit on a slide
  #Publishing here so can verify its right and any changes can be made to both versions
  
  
  #### EM Function ####
  EM_function <- function(X, m, n = input$n){
    if(is.null(n) || missing(n) || is.na(n)) n = 50  #Prevents Crash
    # performs the Expectation-maximization algorithm on data X
    # with m mixing components, and n = 50 iterations at most
    
    X = as.numeric(X); N = length(X)
    
    if (m != 1){
      q=(1:(m-1))*1/m
      quantiles=c(unname(quantile(X,q)), max(X))
      tempmu=numeric(m); tempsigma=numeric(m)
      
      # first part
      tempmu[1]=mean(X[X<=quantiles[1]]); tempsigma[1]=sd(X[X<=quantiles[1]])
      for (i in 2:m){
        tempmu[i]=mean(X[X>quantiles[i-1] & X<=quantiles[i]]); tempsigma[i]=sd(X[X>quantiles[i-1] & X<=quantiles[i]])
      }
    } else {
      tempmu=c(mean(X)); tempsigma=c(sd(X))
    }
    
    # initialize the parameters
    pi = matrix(0,n,m); mu = matrix(0,n,m); sigma = matrix(0,n,m)
    pi[1,] = rep(1/m,m) # weight them equally
    # mu[1,] = seq(from = min(X), to = max(X), length.out = m) # space them out evenly
    # sigma[1,] = rep(diff(range(X))/(6*m),m) # spread them out evenly
    mu[1,] = tempmu
    sigma[1,] = tempsigma
    
    theta = list(pi = pi, mu = mu, sigma = sigma, L = 0)
    
    iter = 1
    while (iter < n){
      # calculate the conditional probabilities
      p = matrix(0,N,m)
      for (j in 1:m){
        p[,j] = theta$pi[iter,j]*dnorm(X, mean = theta$mu[iter,j], sd = theta$sigma[iter,j])
      }
      if (any(is.infinite(p))){
        theta$pi=theta$pi[1:iter,]
        theta$mu=theta$mu[1:iter,]
        theta$sigma=theta$sigma[1:iter,]
        return(theta)
      }
      
      if (iter == 1){
        LogLikelihood = c(sum(log(rowSums(p))))
      } else {
        LogLikelihood = c(LogLikelihood, sum(log(rowSums(p))))
      }
      theta$L = LogLikelihood[length(LogLikelihood)]
      
      p.hat = t(apply(p, 1, function(x){x/sum(x)}))
      if (m == 1){p.hat = t(p.hat)} # R doesn't like mathematicians
      
      # calculate the new parameters
      
      new_pi = apply(p.hat,2,sum)/N
      new_mu = apply(p.hat,2,function(x){sum(x*X)/sum(x)})
      new_sigma = sqrt(apply(p.hat,2,function(x){sum(x*(X-sum(x*X)/sum(x))^2)/sum(x)}))
      
      theta$pi[iter+1,]=new_pi
      theta$mu[iter+1,]=new_mu
      theta$sigma[iter+1,]=new_sigma
      
      if (iter > 1){
        if(abs((LogLikelihood[iter] - LogLikelihood[iter-1])/LogLikelihood[iter-1]) < 1e-5){
          theta$pi=theta$pi[1:(iter+1),]
          theta$mu=theta$mu[1:(iter+1),]
          theta$sigma=theta$sigma[1:(iter+1),]
          break
        }
      }
      
      iter = iter + 1
    }
    
    return(theta)
  }
