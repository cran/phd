doubleres <- function(y, X, X1, nperm=2E4, lambda="lambda.min",flip="FALSE",nfolds=10){
  
  w=nperm
  
  if(!missing(X1)){
    if(length(X1)!=nrow(X)){
      stop('X and X1 should both have n rows')
    }
    
    X <- cbind(X1,X)
  } 
  
  if(length(y)!=nrow(X)){
    stop('the length of y should equal the nr of rows of X')
  }
  

  if(is.numeric(lambda)==TRUE){
    if(lambda<=0){
      stop('lambda should be larger than 0')
    }
  } 
  
  if(!(flip==TRUE|flip==FALSE)){
    stop('argument flip should be TRUE or FALSE')
  } 
  
  if(w<100){
    message('taking the nr of perms w very small leads to inaccurate p-values')
  }
  
  n=nrow(X)
  q = ncol(X)-1
  pen=0     #penalty type
 

  range = 10^seq(5,-5,-0.1)  #range of candidate values for lambda

  
  #STANDARDIZE THE DATA
  X = apply(X, 2, function(x){(x-mean(x))/sd(x)})
  y = (y-mean(y))/sd(y)
  
  
  modelzy <- cv.glmnet(x=X[,-1],y=y,family="gaussian",alpha=pen,intercept=FALSE,lambda=range,nfolds=nfolds) 
  modelzx  <- cv.glmnet(x=X[,-1],y=X[,1],family="gaussian",alpha=pen,intercept=FALSE,lambda=range,nfolds=nfolds)
 
  lam=lambda  #penalty for regression of y on z
  la=lambda   #penalty for regression of x on z
  if(lambda=="lambda.min"){
    lam =  n*modelzy$lambda.min #multiply by n, since glmnet scales the penalty differently
    la = n*modelzx$lambda.min 
  }
  if(lambda=="lambda.1se"){
    lam =  n*modelzy$lambda.1se  
    la = n*modelzx$lambda.1se
  }


  X0 = X[,-1]
  yhat = X0%*%solve( t(X0)%*%X0+lam*diag(q) )%*%t(X0)%*%y      
  xhat = X0%*%solve( t(X0)%*%X0+lam*diag(q) )%*%t(X0)%*%X[,1]  
 
  resids = y-yhat
  resx = X[,1] - xhat
 
 
  T <- numeric(w)
  T[1] <- cor(resx,y)        
 
  for(j in 2:w){
   
    Ystar = yhat + sample(resids,size=n,replace=FALSE)
   
    if(flip==TRUE){
      flips <- 2*rbinom(n,1,0.5)-1
      Ystar <- yhat + flips*resids         #rnorm(n, mean = predictzy, sd=sigma)
    } 
   
    T[j] <- cor(resx,Ystar)        
   
  }
 
  p1 <- sum(T>=T[1])/w 
  p2 <- sum(T<=T[1])/w
  p = 2*min(p1,p2)
 
  p
}
 
   