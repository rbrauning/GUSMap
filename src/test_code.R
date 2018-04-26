setwd("GUSMap/src")
dyn.load("score.so")
library(GUSMap)
data <- simFS(rVec_f = 0.01,epsilon = 0.01,config = c(1,2,1,4,1),meanDepth = 5,nInd = 10,NoDS = 1,seed2 = 79847)

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

tt <- rf_est_FS(depth_Ref = depth_Ref,depth_Alt = depth_Alt,OPGP=OPGP, method="optim")


score_fs_mp_scaled_err <- function(para,depth_Ref,depth_Alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam){
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
    score = score + .Call("score_fs_scaled_err_c",r,para[nSnps],depth_Ref[[fam]],depth_Alt[[fam]],
                          Kaa[[fam]],Kab[[fam]],Kbb[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
  return(score)
}

score_fs_mp_scaled_err(c(GUSMap:::logit2(tt[[1]]),GUSMap:::logit(tt[[2]])),depth_Ref,depth_Alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam=1)


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
                   fn=ll_fs_mp_scaled_err,gr=score_fs_mp_scaled_err,
                   method="BFGS",
                   depth_Ref=depth_Ref,depth_Alt=depth_Alt,bcoef_mat=bcoef_mat,Kab=Kab,
                   nInd=nInd,nSnps=nSnps,OPGP=OPGP,noFam=noFam,seqErr=T)

