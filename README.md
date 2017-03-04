# bayesian-MCMC-modeling
hierarchical bayesian linear regression(HB regression) with gibbs sampling

#load data
yxdata <- read.table("C:/GAUSS/sample/data/hreg_yxdata.txt", quote="\"", comment.char="")
zdata <- read.table("C:/GAUSS/sample/data/hreg_zdata.txt", quote="\"", comment.char="")

y=as.matrix(yxdata[,2])
x=as.matrix(yxdata[,3:ncol(yxdata)])
z=as.matrix(zdata)

n=nrow(y) ## num of data
k=ncol(x) ##column of x data
l=ncol(z) ##column of z data

thdim=k*l ##dimension of theta
ztz=t(z)%*%z

## model:   y_hi = x_hi'beta_i + e_hi;   e_hi ~ N(0, sigma)     */
##          beta_i = theta'z_i + d_i; d_i ~ N(0, delta)         */
##          h: subject or consumer,     i: purchase occasion    */

## data handling
nsubx=function(data){
  data=as.matrix(data)
  n=nrow(data)
  nsub=1				    # num of subjects @
  yrows=as.matrix(1,1,1)				    # num of samples for each subjects @
  id=matrix(data[1],1,1)

   for (i0 in 2:n) {
   i=i0
    if(data[i]==data[i-1]) {yrows[nsub]=yrows[nsub]+1}
    else{
      nsub=nsub+1
      yrows=rbind(yrows,1)
      id=rbind(id,data[i])
    }
   }

  lxy=matrix(0,nsub,1)
  uxy=lxy
    
  for (i0 in 1:nsub) {
    i<-i0
    if(i==1){
      lxy[1]=1
      uxy[1]=yrows[1]
    } else if(i>1) {
      lxy[i]=sum(yrows[1:(i-1)])+1
      uxy[i]=sum(yrows[1:i])
    }
    
  }
  return(list(nsub=nsub,id=id,yrows=yrows,lxy=lxy,uxy=uxy))
}
  result=nsubx(yxdata[,1])
  nsub=result[[1]]
  id=result[[2]]
  yrows=result[[3]]
  lxy=result[[4]]
  uxy=result[[5]]
#  nsub = num of subject (or consumer)
#  id = id number
#  yrows = num of each subject data
#  lxy, uxy = row numbers of first and last data in each subject
#  x[lxy[i]:uxy[i].,] => get x_i (x data of subject i)

#  mcmc iteration  ***/
  nblow=500     # burn in
  smcmc=500
  
#  prior  ***/
    # beta_i ~ N_k(theta'z_i',delta)  @
    
    #  sigma ~ IG(s0,r0) @
  s0=2
  r0=2
  r0n=r0+n
  
  #  theta ~ N_thdim(u0,v0)  @
  u0=matrix(0,thdim,1)
  v0=diag(1,thdim)*100
  v0i=solve(v0)
  v0iu0=v0i%*%u0
  
  #  delta ~ IW_k(f0,g0)  @
  f0=k+2             # degree of freedom
  f0n=f0+nsub
  g0i=diag(1,k)*f0     # scale matrix
  
  #  initial value and box for gibbs sampler  ***/
  beta=matrix(0,nsub,k)         # initial value of all subject
  betam=beta                 # use to evaluate posterior mean
  betas=beta                 # use to evaluate posterior s.d
  betag=matrix(0,smcmc,k)      # use to check convergence
  
  sigma=1                    # sigma is common in all subjects
  sigmai=1/sigma
  sigmag=matrix(0,smcmc,1)
  
  theta=matrix(0,l,k)
  thetag=matrix(0,smcmc,thdim)
  
  library(MCMCpack) # a package provides the random inverse wishart distrbution
  library(matrixcalc) # a package provides matrice transformation, especially transform between vector and matrix
  library(bayesm)
  
  delta=diag(1,k)
  deltai=solve(delta)
  deltag=matrix(0,smcmc,nrow(vech(delta)))
  
# multivariate regression
getmulreg = function(yd,xd,xdtxd,parmat,var,vari,v0i,v0iu0,f0n,g0i){
  rp    = nrow(parmat)
  cp		= ncol(parmat);
  pdim	= rp*cp
  
  vb12	= ( kronecker(xdtxd,vari) + v0i)
  ubn		= kronecker( t(xd),vari )%*%vec(t(yd)) + v0iu0
  par		= solve(vb12)%*%(ubn + t(chol(vb12))%*%rnorm(pdim,0,1))
  parmat	= matrix(par,rp,cp,byrow = TRUE)
  resid 	= yd - xd%*%parmat
  gni		  = g0i + t(resid)%*%resid
  gn  		= solve(gni)
  var    = riwish(f0n,gni)
  vari   = solve(var)
  return(list(parmat=parmat,var=var,vari=vari))
}

# univariate regression
hbreg=function(y=y,x=x,z=z,ztz=ztz,beta=beta,sigma=sigma,sigmai=sigmai,theta=theta,delta=delta,deltai=deltai,s0=s0,s0n=s0n){
  #  beta_i ~ N_k(bm_i,bv_i)  @
  #  sigma ~ IG(sn,r0n)  @
    zt=z%*%theta     # prior of beta => beta ~ N(zt[i,.]',delta)
   
    sse=0            # sum of squared errors (use to generate sigma)
    for (i0 in 1:nsub) {
      i=i0
      yi=y[lxy[i]:uxy[i]]                   # y data of subject i
      xi=x[lxy[i]:uxy[i],]                  # x data of subject i
      bv=solve(deltai+as.numeric(sigmai)*(t(xi)%*%xi))          # conditional posterior variance
      
      bm=bv%*%(deltai%*%(zt[i,])+as.numeric(sigmai)*(t(xi)%*%yi))   # conditional posterior mean
      betai=t(chol(bv))%*%rnorm(k,0,1)+bm            # generate beta_i
      beta[i,]=t(betai)                       # update beta matrix
      
      res=yi-xi%*%betai                        # residuals of subject i
      sse=sse+t(res)%*%res                        # calculate sum of squared errors and sum in all subject
    }
    sn=s0+sse      # posterior scale
    sigma=sn/(rgamma(1,shape=r0n/2,scale=2))     # generate and update sigma
    sigmai=1/sigma                    # update sigmai=1/sigma
    
    #  theta & delta  @
    
    
    result=getmulreg(beta,z,ztz,theta,delta,deltai,v0i,v0iu0,f0n,g0i)
    theta=result$parmat
    delta=result$var
    deltai=result$vari
   
    # input = (data, parameter, prior),  output = (updated parameter) 
    return(list(beta=beta,sigma=sigma,sigmai=sigmai,theta=theta,delta=delta,deltai=deltai))
}

#  do MCMC  ***/
  for (i1 in 1:nblow){
    imcmc=i1
    result=hbreg(y=y,x=x,z=z,ztz=ztz,beta=beta,sigma=sigma,sigmai=sigmai,theta=theta,delta=delta,deltai=deltai,s0=s0,s0n=s0n)
    beta=result$beta
    sigma=result$sigma
    sigmai=result$sigma
    theta=result$theta
    delta=result$delta
    deltai=result$deltai
  } 



for (i1 in 1:smcmc){
  imcmc=i1
  result=hbreg(y=y,x=x,z=z,ztz=ztz,beta=beta,sigma=sigma,sigmai=sigmai,theta=theta,delta=delta,deltai=deltai,s0=s0,s0n=s0n)
  beta=result$beta
  sigma=result$sigma
  sigmai=result$sigma
  theta=result$theta
  delta=result$delta
  deltai=result$deltai
  
  betam=betam+beta               # for posterior mean
  betas=betas+beta^2             # for posterior s.d
  betag[imcmc,]=t(colMeans(beta))
  
  sigmag[imcmc,]=sigma
  
  thetag[imcmc,]=t(vec(t(theta)))   # save mcmc sample
  deltag[imcmc,]=t(vech(delta))
} 

#  posterior  ***/
library(matrixStats)      # Matrice process package
betam=betam/smcmc      # posterior mean
betas=sqrt(abs(betas-smcmc*betam^2)/smcmc);     # posterior s.d

sigmam=colMeans(sigmag);
sigmas=colSds(sigmag);

thetam=matrix(colMeans(thetag),l,k,byrow = TRUE)     # posterior mean
thetas=matrix(colSds(thetag),l,k,byrow = TRUE)      # posterior sd

deltam=xpnd(colMeans(deltag))     # xpnd(w) : create symmetric matrix from w=vech(z)
deltas=xpnd(colSds(deltag))

# convergence check & parameter distribution   
par(mfrow=c(2,4))
plot(thetag[,1],type = "l")
hist(thetag[,1])
plot(thetag[,2],type = "l")
hist(thetag[,2])
plot(deltag[,1],type = "l")
hist(deltag[,1])
plot(deltag[,2],type = "l")
hist(deltag[,2])

par(mfrow=c(2,2))
plot(betag[,1],type = "l")
hist(betag[,1])
plot(betag[,2],type = "l")
hist(betag[,2])

#  output  ***/
print ("===  beta  ===");
print ("mean                 std");
print (cbind(colMeans(betam),colSds(betam)));

print ("===  sigma ===");
print ("mean                 std");
print (cbind(sigmam,sigmas));

print ("===  theta  ===");
print ("mean");
print (thetam);
print ("std");
print (thetas);

print ("===  delta  ===");
print ("mean");
print (deltam);
print ("std");
print (deltas);
