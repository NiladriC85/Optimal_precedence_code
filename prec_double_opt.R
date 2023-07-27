rm(list=ls())

library(cubature)
library(foreach)
library(iterators)
library(GA)
library(knitr)
library(tidyverse)
library(hrbrthemes)

  m = 75
  n = 64
  gamma = 4
  r=2
alpha = 0.05
eta = 2
prec.prob = function(k){
  result1 = choose(k+r-1,r-1)*choose(m+n-k-r,n-r)/choose(m+n,m)
  return(result1)
}

############################################### obtaining n1 value #######################################

f1 = function(x) {sum(prec.prob(c(x:m)))}
vec = c(1:m)
a = min(which(sapply(vec, f1)<alpha)) ### this is same as c
c=a

f2 = function(n1) sum(sapply(c(a:m), function(k) choose(k+r-1,r-1)*choose(m+n1-k-r,n1-r)/choose(m+n1,m)))


vec1 = c(1:n)
n1.1 = max(which(sapply(vec1, f2)> eta*f1(a) )) ### n1 value for which a = c 
n1.1

f2(n1.1)
## we need to take n1 greater than n1.1

######## obtaining the b value and the power #####################################

n1 = n1.1+1 ## initial value for n1

################################################ alpha value ##########################################
prob0 = function(b){
  sum1 = 0
  for(i in 1:(n-n1+1)){
    d = i-1
    for(k1 in seq((r+d):n)){
      for(k in seq(r:k1)){
        sum1 = sum1+(choose(b-c-1+k1-k,k1-k)*choose(k+c-1,k)*choose(m+n-b-k1,m-b)/choose(m+n,m))*(choose(r+d-1,d)*choose(n-r-d,n-n1-d)/choose(n,n1))
      }
    }
  }
  result11 = 1-sum1
  return(result11)
}

################################################ power calculation ##################################
b = c+1 ### initial assignment of b

### first we calculate type II probability for nonzero d ############
type2.prob = function(d){
  prob1 = function(u1){
    fun = function(x1){
      (r+d)*choose(n,r+d)*u1^(gamma*(r+d))*x1^(r+d-1)*(1-x1*u1^gamma)^(n-r-d)
    }
    result1 = pcubature(fun,lowerLimit = c(0),upperLimit = c(1),tol = 1e-05)$integral
    return(result1)
  }
  
  prob2 = function(z4,u1){
    fun = function(x){
      A = r*d*choose(r+d,r)*choose(n,r+d)*u1^(gamma*(r+d))*(z4^gamma+x[1]*x[2]*(1-z4^gamma))^(r-1)
      B = (x[2])^d*(1-z4^gamma)^(d-1)*(1-x[1])^(d-1)*(1-z4^gamma*u1^gamma-x[2]*u1^gamma*(1-z4^gamma))^(n-r-d)*(1-z4^gamma)^2
      result2 = A*B
      return(result2)
    }
    result3 = pcubature(fun,lowerLimit = c(0,0),upperLimit = c(1,1),tol = 1e-05)$integral
    return(result3)
  }
  ## we define the probability functions below
  prob11 = function(u1){
    result4 = prob1(u1)*b*choose(m,b)*u1^(b-1)*(1-u1)^(m-b)
    return(result4)
  } 
  
  prob21 = function(x2){
    result5 = prob2(x2[1],x2[2])*choose(m,b)*(b-c+1)*(b-c)*choose(b,c-1)*x2[1]^(c-1)*x2[2]^(b-1)*(1-x2[1])^(b-c-1)*(1-x2[2])^(m-b)
    return(result5)
  }
  #prob2(0.9,0.9)
  
  A1 = pcubature(prob11,lowerLimit = c(0),upperLimit = c(1),tol = 1e-05)$integral
  B1 = pcubature(prob21,lowerLimit = c(0,0),upperLimit = c(1,1),tol = 1e-05)$integral
  Result6 = A1-B1
  return(Result6)
}

vec11 = c(1:(n-n1))

#type2.prob((n-n1))

#sum(sapply(vec11, type2.prob))

### then we calculate type II probability for d = 0 ############

prob3 = function(x){
  Result7 = r*choose(n,r)*c*choose(m,c)*(x[1]*x[2]^gamma)^(r-1)*(1-x[1]*x[2]^gamma)^(n-r)*x[2]^(gamma+c-1)*(1-x[2])^(m-c)
  return(Result7)
}

#pcubature(prob3,lowerLimit = c(0,0),upperLimit = c(1,1),tol = 1e-05)$integral

### now we will define the objective function ###################################################
### note that n1 > n1.1

power = function(n1){
  #n=85
  ## first we calculate the b value using the alpha function define above
  m = 75
  gamma = 4
  r=4
  n=61
  vec2 = c((c+1):m)
  unit1 = min(which(sapply(vec2, prob0)<alpha))
  b = c+unit1
  # prob0(b) ### now we calculate the power
  vec11 = c(1:(n-n1))
  pwr = 1-(pcubature(prob3,lowerLimit = c(0,0),upperLimit = c(1,1),tol = 1e-05)$integral+sum(sapply(vec11, type2.prob)) )
  return(pwr)
}

#power(n1.1+10)

objective.fn = function(n1){
  n12 = n1+2
  result8 = 1-(abs(power(n12)-power(n1)))^(ifelse(f1(a)<alpha*(1+epsilon) & f1(a)>alpha*(1-epsilon),1,0) )
  return(result8)
}

#objective.fn(81)

vec3 = c((n1.1+1):(n-2))

min(which.max(sapply(vec3, objective.fn)))
power(n1.1)
n1.1

lbound = n1.1
ubound = n
GA1 <- ga(type = "real-valued", fitness = objective.fn, 
          lower = lbound, upper = ubound,
          popSize = (n-n1.1+1), maxiter = 1000, run = 100)


summary(GA1)

plot(GA1)

