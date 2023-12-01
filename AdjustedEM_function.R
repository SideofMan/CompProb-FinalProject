  #I wanted to show the EM funtion in our presentation but the original version was a little lengthy, this does the exact same thing just a little more condensed 
  #So I can fit on a slide
  #Publishing here so can verify its right and any changes can be made to both versions
  
  
  #### EM Function ####
  EM_function <- function(X, m, n = input$n){ #performs the Expectation-maximization algorithm on data X
    if(is.null(n) || missing(n) || is.na(n)) n = 50   #Prevents Crash
    # with m mixing components, and n = 50 iterations at most
    X = as.numeric(X); N = length(X)

    #Initialize the parameters
    pi = matrix(0,n,m); mu = matrix(0,n,m); sigma = matrix(0,n,m)
    pi[1,] = rep(1/m,m) #Weight them equally
    if (m != 1){ #IF m>1
      bins = cut(X, breaks = c(-Inf,quantile(X,(1:m)/m)), include.lowest = TRUE, labels = FALSE) #Generate m quantiles
      mu[1,] = tapply(X, bins, mean); sigma[1,] = tapply(X, bins, sd) #Init mu,sigma at mean/sd of each quantile
    } else {
      mu[1,]=c(mean(X)); sigma[1,]=c(sd(X))
    }
    theta = list(pi = pi, mu = mu, sigma = sigma)
    
    iter = 1
    while (iter < n){
      # calculate the conditional probabilities
      p = sapply(1:m, function(j) theta$pi[iter,j] * dnorm(X, mean = theta$mu[iter,j], sd = theta$sigma[iter,j]))
      
      if (any(is.infinite(p))) return(lapply(theta, function(x) x[1:iter, ])) #End condition
      
      p.hat = matrix(t(apply(p, 1, function(x) x / sum(x))), ncol = m)
      
      #Calculate the new parameters
      theta$pi[iter+1,] = apply(p.hat,2,sum)/N
      theta$mu[iter+1,] = apply(p.hat,2,function(x){sum(x*X)/sum(x)})
      theta$sigma[iter+1,]= sqrt(apply(p.hat,2,function(x){sum(x*(X-sum(x*X)/sum(x))^2)/sum(x)}))
      
      #Calculate Log Likelihood
      if (iter == 1){#Init
        LogLikelihood = c(sum(log(rowSums(p))))
      } else { 
        LogLikelihood = c(LogLikelihood, sum(log(rowSums(p))))
        if(abs((LogLikelihood[iter] - LogLikelihood[iter-1])/LogLikelihood[iter-1]) < 1e-5){ #End condition
          theta <- lapply(theta, function(x) x[1:(iter + 1), ])
          break
        }
      }
      iter = iter + 1
    }
    
    return(theta)
  }
