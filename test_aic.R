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

df = read.csv(file = "sim-data.csv")
x = df[,1]
modes = 1:10
a = c()
b = c()
n = 50
for(m in modes) {
  EM_fit = EM_function(x, m)
  theta = EM_fit
  if (m == 1){
    theta$pi = matrix(theta$pi, ncol = 1)
    theta$mu = matrix(theta$mu, ncol = 1)
    theta$sigma = matrix(theta$sigma, ncol = 1)
  }
  final=nrow(theta$pi)
  
  
  x_values = seq(min(x),max(x),by=0.05)
  y_values = matrix(0, length(x_values), m)
  
  # initial guesstimate
  for (j in 1:m){
    y_values[,j] = theta$pi[1,j]*dnorm(x_values, mean = theta$mu[1,j], sd = theta$sigma[1,j])
  }
  y_values_initial = rowSums(y_values)
  
  # output at user defined step
  EM_step=m
  for (j in 1:m){
    y_values[,j] = theta$pi[EM_step,j]*dnorm(x_values, mean = theta$mu[EM_step,j], sd = theta$sigma[EM_step,j])
  }
  y_values_EM_step = rowSums(y_values)
  
  # final output
  final=nrow(theta$pi)
  for (j in 1:m){
    y_values[,j] = theta$pi[final,j]*dnorm(x_values, mean = theta$mu[final,j], sd = theta$sigma[final,j])
  }
  y_values = rowSums(y_values)
  
  #calculate AIC 
  summedValues = sum(log(y_values))
  a = c(a, 2 * m - 2 * summedValues)
  b = c(b, m * log(n) - 2 * summedValues)
}

plot(b, type = 'b', col = "red", xlab = "Modes", ylab = "AIC/BIC")
points(a, type = 'b', col = "blue")
legend(8,840,legend=c("AIC","BIC"), fill = c("blue","red") )
