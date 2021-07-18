library(bbmle)

source(paste('path-to/primary-functions-1.R',sep=''))

### DATA ENTRY ########################################

b = # bid data goes here
sold = # auction outcome data goes here

auctions = length(b)

min.b = min(b)
max.b = max(b)
seq.b = seq(min.b,max.b,len=101)

### STAGE 1 ESTIMATION ########################################

k1 = 2

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
  print(p)
  
  d = integrate(fhat0.stage.one,unname(p['min.x']),unname(p['max.x']),p=p,k=k1)$value
  
  -sum(log(v.fhat.stage.one(b,p,k1,d)))/auctions
}
parnames(NLL.stage.one) = names(starting.guess)
fit.stage.one = mle2(NLL.stage.one,start=starting.guess,method='L-BFGS-B','optim',lower=constraints.min,upper=constraints.max,control=list(maxit=10000))
summary(fit.stage.one)
params.stage.one = coef(fit.stage.one)
d.stage.one = integrate(fhat0.stage.one,unname(params.stage.one['min.x']),unname(params.stage.one['max.x']),p=params.stage.one,k=k1)$value

plot(seq.b,sapply(seq.b,function(i) Fhat.stage.one(i,params.stage.one,k1,d.stage.one)),ylim=c(0,1))
lines(ecdf(b))

### STAGE 2 ESTIMATION ########################################

y = b[sold==1]
x = b[sold==0]

min.y = min(y)
max.x = max(x)

k2 = 2

min.std = 0.1
max.std = 1e6
min.mean = -1e6
max.mean = 1e6
min.params = -1e6
max.params = 1e6
constraints.min = eval(parse(text=lower.stage.two(k2,min.params,min.mean,min.std)))
constraints.max = eval(parse(text=upper.stage.two(k2,max.params,max.mean,max.std)))
starting.guess = eval(parse(text=starter.stage.two(k2,mean(c(mean(x),mean(y))),mean(c(sd(x),sd(y))))))
NLL.stage.two = function(p) {
  print(p)
  
  d = integrate(fhat0.stage.two,0,Inf,p=p,k=k2)$value
  integral = fhat.G.integral.stage.two(params.stage.one,p,k1,k2,d.stage.one,d)
  
  -sum(log(v.f.bids.greater.stage.two(y,params.stage.one,p,k1,k2,d.stage.one,d,integral)))-sum(log(v.f.bids.lesser.stage.two(x,params.stage.one,p,k1,k2,d.stage.one,d,integral)))
}
parnames(NLL.stage.two) = names(starting.guess)
fit.stage.two = mle2(NLL.stage.two,start=starting.guess,method='L-BFGS-B','optim',lower=constraints.min,upper=constraints.max,control=list(maxit=10000))
summary(fit.stage.two)
params.stage.two = coef(fit.stage.two)
d.stage.two = integrate(fhat0.stage.two,0,Inf,p=params.stage.two,k=k2)$value
integral.stage.two = fhat.G.integral.stage.two(params.stage.one,params.stage.two,k1,k2,d.stage.one,d.stage.two)

plot(seq.b,sapply(seq.b,function(i) Fhat.stage.two(i,params.stage.two,k2,d.stage.two)),ylim=c(0,1))

