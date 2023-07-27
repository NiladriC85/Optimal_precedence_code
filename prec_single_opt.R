library(cubature)
library(foreach)
library(iterators)
library(GA)
library(knitr)
library(tidyverse)
library(hrbrthemes)

#gamma = 1
alpha = 0.05

power1 = function(m,a,n,r,gamma){
  fun <- function(x) {
    A = (x[2]^gamma+x[1]*(1-x[2]^gamma))^(r-1)*((1-x[2]^gamma)*(1-x[1]))^(n-r)
    B = x[2]^(a-1)*(1-x[2])^(m-a)*(1-x[2]^gamma)
    m*n*choose((n-1),(r-1))*choose((m-1),(a-1))*A*B
  }
  
  #result1 = adaptIntegrate(fun,lowerLimit = c(0,0),upperLimit = c(1,1))$integral
  result1 = pcubature(fun,lowerLimit = c(0,0),upperLimit = c(1,1),tol = 1e-05)$integral
  return(result1)
}

##########################3 plotting the power function ##########################################
objective.fn1 = function(n){
  m =100
  gamma = 2
  r = 2 
  prob = function(k){
    result11 = choose(k+r-1,k)*choose(m-k+n-r,n-r)/choose(m+n,n)
    return(result11)
  } 
  f1 = function(x) {sum(prob(c(x:m)))}
  #  f1(20)
  vec = c(1:m)
  a = min(which(sapply(vec, f1)<alpha))
  #  f1(a)
  result2 = power1(m,a,n,r,gamma)
  return(result2)
}

objective.fn2 = function(n){
  m =100
  gamma = 2
  r = 3 
  prob = function(k){
    result11 = choose(k+r-1,k)*choose(m-k+n-r,n-r)/choose(m+n,n)
    return(result11)
  } 
  f1 = function(x) {sum(prob(c(x:m)))}
  #  f1(20)
  vec = c(1:m)
  a = min(which(sapply(vec, f1)<alpha))
  #  f1(a)
  result2 = power1(m,a,n,r,gamma)
  return(result2)
}

objective.fn3 = function(n){
  m =100
  gamma = 2
  r = 4 
  prob = function(k){
    result11 = choose(k+r-1,k)*choose(m-k+n-r,n-r)/choose(m+n,n)
    return(result11)
  } 
  f1 = function(x) {sum(prob(c(x:m)))}
  #  f1(20)
  vec = c(1:m)
  a = min(which(sapply(vec, f1)<alpha))
  #  f1(a)
  result2 = power1(m,a,n,r,gamma)
  return(result2)
}

xlim = c(10:99)
yvalue1 = sapply(c(10:99),objective.fn1)
yvalue2 = sapply(c(10:99),objective.fn2)
yvalue3 = sapply(c(10:99),objective.fn3)

data = data.frame(xlim,ylim=c(yvalue1,yvalue2,yvalue3),group=rep(c("r=2","r=3","r=4"),each=90))

ggplot(data,aes(x=xlim,y=ylim,linetype=group))+geom_line(color="black")+theme_classic()+xlab("n")+ylab("power")

##############################################################################################
epsilon = 0.2
objective.fn = function(n){
#n = 130
  m =100
  gamma = 3
  r = 2
  prob = function(k){
    result11 = choose(k+r-1,k)*choose(m-k+n-r,n-r)/choose(m+n,n)
    return(result11)
  } 
  f1 = function(x) {sum(prob(c(x:m)))}
#  f1(20)
  vec = c(1:m)
  a = min(which(sapply(vec, f1)<alpha))
#  f1(a)
  n1 = n+2
  result2 = 1-(abs(power1(m,a,n1,r,gamma)-power1(m,a,n,r,gamma)))^(ifelse(f1(a)<alpha*(1+epsilon) & f1(a)>alpha*(1-epsilon),1,0) )
  return(result2)
}



lbound = 40
ubound = 140
GA2 <- ga(type = "real-valued", fitness = objective.fn, 
          lower = lbound, upper = ubound,
          popSize = 80, maxiter = 1000, run = 100)


summary(GA2)

plot(GA2)








