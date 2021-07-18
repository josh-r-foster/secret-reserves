library(bbmle)
library(devtools)

source_url('https://github.com/josh-r-foster/secret-reserves/raw/main/r-reservation-functions-1.R')

experiments = 500
auctions = 100
k1 = 2
k2 = 2

reserve.prices = matrix(0,ncol=auctions,nrow=experiments)
winning.bids = matrix(0,ncol=auctions,nrow=experiments)

stage.one.parameters = matrix(0,ncol=(k1+4),nrow=experiments)
stage.one.denoms = rep(0,experiments)
stage.two.parameters = matrix(0,ncol=(k2+2),nrow=experiments)
stage.two.denoms = rep(0,experiments)

experiment = 1
while(experiment<=experiments) {
  
  print(paste('Estimating experiment ',as.character(experiment),' of ',as.character(experiments),'...',sep=''))
  
  ### DATA GENERATION ###
  
  # exponential distributions 
  # r = rexp(auctions,1/3)
  b = rexp(auctions,1/3)
  
  # normal distributions
  # r = rnorm(auctions,3,1)
  # b = rnorm(auctions,3,1) 
  
  # uniform distributions 
  r = runif(auctions,2,4)
  # b = runif(auctions,2,4)
  
  min.r = min(r)
  max.r = max(r)
  seq.r = seq(min.r,max.r,len=101)
  
  min.b = min(b)
  max.b = max(b)
  seq.b = seq(min.b,max.b,len=101)
  
  ### STAGE 1 ESTIMATION ###
  
  min.std = 0.1
  max.std = 1e6
  min.mean = -1e6
  max.mean = 1e6
  min.params = -1e6
  max.params = 1e6
  constraints.min = eval(parse(text=lower.stage.one(k1,min.params,min.mean,min.std,max.b)))
  constraints.max = eval(parse(text=upper.stage.one(k1,max.params,max.mean,max.std,min.b)))
  starting.guess = eval(parse(text=starter.stage.one(k1,mean(b),sd(b),min.b,max.b)))
  NLL.stage.one = function(p) {
    d = integrate(fhat0.stage.one,unname(p['min.x']),unname(p['max.x']),p=p,k=k1)$value
    -sum(log(v.fhat.stage.one(b,p,k1,d)))/auctions
  }
  parnames(NLL.stage.one) = names(starting.guess)
  fit.stage.one = tryCatch(
    mle2(NLL.stage.one,start=starting.guess,method='L-BFGS-B','optim',lower=constraints.min,upper=constraints.max,control=list(maxit=10000)),
    error=function(e) e
  )
  if(!inherits(fit.stage.one,'error')){
    print(summary(fit.stage.one))
    params.stage.one = coef(fit.stage.one)
    d.stage.one = integrate(fhat0.stage.one,unname(params.stage.one['min.x']),unname(params.stage.one['max.x']),p=params.stage.one,k=k1)$value
    
    plot(seq.b,sapply(seq.b,function(i) Fhat.stage.one(i,params.stage.one,k1,d.stage.one)),ylim=c(0,1))
    lines(ecdf(b))
  } else {
    print(paste('Experiment ',as.character(experiment),' failed in Stage 1.  Starting over...',sep=''))
    next 
  } 
  
  ### STAGE 2 ESTIMATION ### 
  
  sold = as.numeric(b>=r)
  
  y = b[sold==1]
  x = b[sold==0]
  
  min.y = min(y)
  max.x = max(x)
  
  min.std = 0.2
  max.std = 1e6
  min.mean = -1e6
  max.mean = 1e6
  min.params = -1e6
  max.params = 1e6
  constraints.min = eval(parse(text=lower.stage.two(k2,min.params,min.mean,min.std)))
  constraints.max = eval(parse(text=upper.stage.two(k2,max.params,max.mean,max.std)))
  starting.guess = eval(parse(text=starter.stage.two(k2,mean(c(mean(x),mean(y))),mean(c(sd(x),sd(y))))))
  NLL.stage.two = function(p) {
    d = integrate(fhat0.stage.two,0,Inf,p=p,k=k2)$value
    integral = fhat.G.integral.stage.two(params.stage.one,p,k1,k2,d.stage.one,d)
    -sum(log(v.f.bids.greater.stage.two(y,params.stage.one,p,k1,k2,d.stage.one,d,integral)))-sum(log(v.f.bids.lesser.stage.two(x,params.stage.one,p,k1,k2,d.stage.one,d,integral)))
  }
  parnames(NLL.stage.two) = names(starting.guess)
  
  fit.stage.two = tryCatch(
    mle2(NLL.stage.two,start=starting.guess,method='L-BFGS-B','optim',lower=constraints.min,upper=constraints.max,control=list(maxit=10000)),
    error=function(e) e
  )
  if(!inherits(fit.stage.two,'error')){
    print(summary(fit.stage.two))
    params.stage.two = coef(fit.stage.two)
    d.stage.two = integrate(fhat0.stage.two,0,Inf,p=params.stage.two,k=k2)$value
    integral.stage.two = fhat.G.integral.stage.two(params.stage.one,params.stage.two,k1,k2,d.stage.one,d.stage.two)
    
    plot(seq.r,sapply(seq.r,function(i) Fhat.stage.two(i,params.stage.two,k2,d.stage.two)),ylim=c(0,1))
    lines(ecdf(r))
  } else {
    print(paste('Experiment ',as.character(experiment),' failed in Stage 2.  Starting over...',sep=''))
    next 
  } 
  
  reserve.prices[experiment,] = r
  winning.bids[experiment,] = b
  
  stage.one.parameters[experiment,] = params.stage.one
  stage.one.denoms[experiment] = d.stage.one
  stage.two.parameters[experiment,] = params.stage.two
  stage.two.denoms[experiment] = d.stage.two
  
  experiment = experiment+1
}