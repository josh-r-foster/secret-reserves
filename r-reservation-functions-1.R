### STAGE 1 ESTIMATION ########################################

### PARAMETER SETTING FUNCTIONS ### 

starter.stage.one = function(k,mean,std,min,max) {
  u = 'c('
  for(i in seq(k-1)) {
    u = paste(u,'a',as.character(i),'=0,',sep='')
  }
  u = paste(u,'a',as.character(k),'=0,m=',as.character(mean),',s=',as.character(std),',min.x=',as.character(min),',max.x=',as.character(max),')',sep='')
  return(u)
}

lower.stage.one = function(k,as,m,s,max.b) {
  u = 'c('
  for(i in seq(k-1)) {
    u = paste(u,'a',as.character(i),'=',as.character(as),',',sep='')
  }
  u = paste(u,'a',as.character(k),'=',as.character(as),',m=',as.character(m),',s=',as.character(s),',min.x=0,max.x=',as.character(max.b),')',sep='')
  return(u)
}

upper.stage.one = function(k,as,m,s,min.b) {
  u = 'c('
  for(i in seq(k-1)) {
    u = paste(u,'a',as.character(i),'=',as.character(as),',',sep='')
  }
  u = paste(u,'a',as.character(k),'=',as.character(as),',m=',as.character(m),',s=',as.character(s),',min.x=',as.character(min.b),',max.x=Inf)',sep='')
  return(u)
}

### CORE DENSITY ESTIMATION FUNCTIONS ### 

fhat0.hermite.stage.one = function(k) {
  u = '1'
  for(i in seq(k)) {
    u = paste(u,'+a',as.character(i),'*((x-m)/s)^',as.character(i),sep='')
  }
  u = paste('(',u,')^2*dnorm(x,m,s)/(pnorm(max.x,m,s)-pnorm(min.x,m,s))',sep='')
  return(u)
}

fhat0.stage.one = function(x,p,k) {
  for(i in seq(k)) {
    eval(parse(text=paste('a',as.character(i),' = p[[\'a',as.character(i),'\']]',sep='')))
  }
  m = p[['m']]
  s = p[['s']]
  min.x = p[['min.x']]
  max.x = p[['max.x']]
  eval(parse(text=fhat0.hermite.stage.one(k)))
}

fhat.stage.one = function(x,p,k,d) {
  if(x<unname(p['min.x']-1e-7)) {
    return(0)
  } else if(x>unname(p['max.x']+1e-7)) {
    return(0)
  } else {
    return(fhat0.stage.one(x,p,k)/d)
  }
}
v.fhat.stage.one = Vectorize(fhat.stage.one,'x')

Fhat.stage.one = function(x,p,k,d) {
  if(x<unname(p['min.x'])) {
    return(0)
  } else if(x>unname(p['max.x'])) {
    return(1)
  } else {
    return(integrate(v.fhat.stage.one,unname(p['min.x']),x,p=p,k=k,d=d)$value)
  }
}

### STAGE 2 ESTIMATION ########################################

### PARAMETER SETTING FUNCTIONS ### 

starter.stage.two = function(k,mean,std) {
  u = 'c('
  for(i in seq(k-1)) {
    u = paste(u,'a',as.character(i),'=0,',sep='')
  }
  u = paste(u,'a',as.character(k),'=0,m=',as.character(mean),',s=',as.character(std),')',sep='')
  return(u)
}

lower.stage.two = function(k,as,m,s) {
  u = 'c('
  for(i in seq(k-1)) {
    u = paste(u,'a',as.character(i),'=',as.character(as),',',sep='')
  }
  u = paste(u,'a',as.character(k),'=',as.character(as),',m=',as.character(m),',s=',as.character(s),')',sep='')
  return(u)
}

upper.stage.two = function(k,as,m,s) {
  u = 'c('
  for(i in seq(k-1)) {
    u = paste(u,'a',as.character(i),'=',as.character(as),',',sep='')
  }
  u = paste(u,'a',as.character(k),'=',as.character(as),',m=',as.character(m),',s=',as.character(s),')',sep='')
  return(u)
}

### CORE DENSITY ESTIMATION FUNCTIONS ### 

fhat0.hermite.stage.two = function(k) {
  u = '1'
  for(i in seq(k)) {
    u = paste(u,'+a',as.character(i),'*((x-m)/s)^',as.character(i),sep='')
  }
  u = paste('(',u,')^2*dnorm(x,m,s)/(1-pnorm(0,m,s))',sep='')
  return(u)
}

fhat0.stage.two = function(x,p,k) {
  for(i in seq(k)) {
    eval(parse(text=paste('a',as.character(i),' = p[[\'a',as.character(i),'\']]',sep='')))
  }
  m = p[['m']]
  s = p[['s']]
  eval(parse(text=fhat0.hermite.stage.two(k)))
}

fhat.stage.two = function(x,p,k,d) {
  fhat0.stage.two(x,p,k)/d
}
v.fhat.stage.two = Vectorize(fhat.stage.two,'x')

Fhat.stage.two = function(x,p,k,d) {
  integrate(v.fhat.stage.two,0,x,p=p,k=k,d=d)$value
}

### SECRET RESERVATION DISTRIBUTION FUNCTIONS ### 

fhat.G.stage.two = function(y,p1,p2,k1,k2,d1,d2) { 
  sapply(y,function(z) {fhat.stage.one(z,p1,k1,d1)*Fhat.stage.two(z,p2,k2,d2)}) 
}
fhat.G.integral.stage.two = function(p1,p2,k1,k2,d1,d2) {
  integrate(fhat.G.stage.two,unname(params.stage.one['min.x']),unname(params.stage.one['max.x']),p1=p1,p2=p2,k1=k1,k2=k2,d1=d1,d2=d2)$value
}

f.bids.greater.stage.two = function(x,p1,p2,k1,k2,d1,d2,integral) {
  fhat.stage.one(x,p1,k1,d1)*Fhat.stage.two(x,p2,k2,d2)/integral
}
v.f.bids.greater.stage.two = Vectorize(f.bids.greater.stage.two,'x')

f.bids.lesser.stage.two = function(x,p1,p2,k1,k2,d1,d2,integral) {
  fhat.stage.one(x,p1,k1,d1)*(1-Fhat.stage.two(x,p2,k2,d2))/(1-integral)
}
v.f.bids.lesser.stage.two = Vectorize(f.bids.lesser.stage.two,'x')