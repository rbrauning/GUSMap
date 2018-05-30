setwd("GUSMap/src")
dyn.unload("score.so")
dyn.load("score.so")
library(GUSMap)

config <- sample(c(1,2,4),size = 10,replace = T, prob = c(1,2,2)/5)

data <- simFS(rVec_f = 0.01,epsilon = 0.01,config = config,meanDepth = 5,nInd = 100,NoDS = 1,seed2 = 10485)

depth_Ref <- list(data$depth_Ref)
depth_Alt <- list(data$depth_Alt)
nInd <- list(as.integer(data$nInd))
nSnps = data$nSnps
OPGP <- list(as.integer(data$OPGP))
noFam = 1

bcoef_mat <- Kab <- vector(mode="list", length=noFam)
for(fam in 1:noFam){
  bcoef_mat[[fam]] <- choose(depth_Ref[[fam]]+depth_Alt[[fam]],depth_Ref[[fam]])
  Kab[[fam]] <- bcoef_mat[[fam]]*(1/2)^(depth_Ref[[fam]]+depth_Alt[[fam]])
}

system.time({tt <- rf_est_FS(depth_Ref = depth_Ref,depth_Alt = depth_Alt,OPGP=OPGP, method="optim",reltol=1e-25)})


score_fs_mp_scaled_err <- function(para,depth_Ref,depth_Alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam,seqErr=T){
  ## untransform the parameters
  r <- GUSMap:::inv.logit2(para[1:(nSnps-1)])
  epsilon = GUSMap:::inv.logit(para[nSnps])
  ## define likelihood
  score = numeric(nSnps)
  # define the density values for the emission probs
  Kaa <- Kbb <- vector(mode = "list", length=noFam)
  for(fam in 1:noFam){
    Kaa[[fam]] <- bcoef_mat[[fam]]*(1-epsilon)^depth_Ref[[fam]]*epsilon^depth_Alt[[fam]]
    Kbb[[fam]] <- bcoef_mat[[fam]]*(1-epsilon)^depth_Alt[[fam]]*epsilon^depth_Ref[[fam]]
  }
  for(fam in 1:noFam)
    score = score + .Call("score_fs_scaled_err_c",r,epsilon,
                          depth_Ref[[fam]],depth_Alt[[fam]], Kaa[[fam]],Kab[[fam]],Kbb[[fam]],
                          OPGP[[fam]],nInd[[fam]],nSnps)
  return(-score)
}
score_fs_mp_scaled_err <- function(para,depth_Ref,depth_Alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam,seqErr=T){
    ## untransform the parameters
    r <- GUSMap:::inv.logit2(para[1:(nSnps-1)])
    if(seqErr)
      epsilon <- GUSMap:::inv.logit(para[nSnps])
    ## define likelihood
    score = numeric(nSnps)
    # define the density values for the emission probs
    Kaa <- Kbb <- vector(mode = "list", length=noFam)
    for(fam in 1:noFam){
      Kaa[[fam]] <- bcoef_mat[[fam]]*(1-epsilon)^depth_Ref[[fam]]*epsilon^depth_Alt[[fam]]
      Kbb[[fam]] <- bcoef_mat[[fam]]*(1-epsilon)^depth_Alt[[fam]]*epsilon^depth_Ref[[fam]]
    }
    if(seqErr){
      score = numeric(nSnps)
      for(fam in 1:noFam)
        score = score + .Call("score_fs_scaled_err_c",r,epsilon,depth_Ref[[fam]],depth_Alt[[fam]],
                            Kaa[[fam]],Kab[[fam]],Kbb[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
    } else{
      score = numeric(nSnps - 1)
      for(fam in 1:noFam)
        score = score + .Call("score_fs_scaled_c",r,Kaa[[fam]],Kab[[fam]],Kbb[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
    }
    return(-score)
}


test <- score_fs_mp_scaled_err(c(GUSMap:::logit2(tt[[1]]),GUSMap:::logit(tt[[2]])),depth_Ref,depth_Alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam=1)

system.time({optim.MLE <- optim(c(GUSMap:::logit2(rep(0.01,data$nSnps-1))),
                   fn=ll_fs_mp_scaled_err,gr=score_fs_mp_scaled_err,
                   method="BFGS",
                   depth_Ref=depth_Ref,depth_Alt=depth_Alt,bcoef_mat=bcoef_mat,Kab=Kab,
                   nInd=nInd,nSnps=nSnps,OPGP=OPGP,noFam=noFam,seqErr=F)})

test2 <- grad( function(u) GUSMap:::ll_fs_mp_scaled_err(u, depth_Ref,depth_Alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam=1,T),
      c(GUSMap:::logit2(tt$rf),GUSMap:::logit(tt$epsilon)))


para <- c(rep(GUSMap:::logit2(0.01),4),GUSMap:::logit(0.01))

GUSMap:::ll_fs_mp_scaled_err(para,depth_Ref,depth_Alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam=1,T)

r <- GUSMap:::inv.logit2(para[1:(nSnps-1)])
epsilon = GUSMap:::inv.logit(para[nSnps])
## define likelihood
score = numeric(nSnps)
# define the density values for the emission probs
Kaa <- Kbb <- vector(mode = "list", length=noFam)
for(fam in 1:noFam){
  Kaa[[fam]] <- bcoef_mat[[fam]]*(1-epsilon)^depth_Ref[[fam]]*epsilon^depth_Alt[[fam]]
  Kbb[[fam]] <- bcoef_mat[[fam]]*(1-epsilon)^depth_Alt[[fam]]*epsilon^depth_Ref[[fam]]
}
.Call("score_fs_scaled_err_c",tt[[1]],tt[[2]],depth_Ref[[fam]],depth_Alt[[fam]],
      Kaa[[fam]],Kab[[fam]],Kbb[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)

ps <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% 1:8)))))[-1] - 1
ms <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% c(1:4,9:12))))))[-1] - 1
npar <- c(length(ps),length(ms))

optim.MLE <- optim(para,ll_fs_mp_scaled_err,method="BFGS",control=optim.arg,
                   depth_Ref=depth_Ref,depth_Alt=depth_Alt,bcoef_mat=bcoef_mat,Kab=Kab,
                   nInd=nInd,nSnps=nSnps,OPGP=OPGP,noFam=noFam,
                   seqErr=seqErr)


optim.MLE <- optim(c(GUSMap:::logit2(rep(0.01,4)),GUSMap:::logit(0.01)),
                   fn=GUSMap:::ll_fs_mp_scaled_err,gr=score_fs_mp_scaled_err,
                   method="BFGS",
                   depth_Ref=depth_Ref,depth_Alt=depth_Alt,bcoef_mat=bcoef_mat,Kab=Kab,
                   nInd=nInd,nSnps=nSnps,OPGP=OPGP,noFam=noFam,seqErr=T)

