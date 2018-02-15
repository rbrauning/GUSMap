##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2018 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################
### Wrapper functions required for calling C code.
### Author: Timothy Bilton
### Date: 06/02/18

#### Wrapper function for likelihoods of the GBS HMM when
#### Likelihood function is written in C in the file 'likelihoods.c'

#' @useDynLib GUSMap
## r.f.'s are equal
ll_fs_mp_scaled_err <- function(para,depth_Ref,depth_Alt,OPGP,nInd,nSnps,noFam,seqErr,samPara){
  ## untransform the parameters
  r <- inv.logit2(para[1:(nSnps-1)])
  if(seqErr)
    epsilon = inv.logit(para[nSnps])
  else
    epsilon = 0
  if(samPara){
    if(epsilon)
      alpha = exp(para[nSnps+1])
    else
      alpha = exp(para[nSnps])
  }
      
  ## define likelihood
  llval = 0
  # define the density values for the emission probs
  bcoef_mat <- Kab <- Kaa <- Kbb <- vector(mode = "list", length=noFam)
  for(fam in 1:noFam){
    bcoef_mat[[fam]] <- choose(depth_Ref[[fam]]+depth_Alt[[fam]],depth_Ref[[fam]])
    if(samPara)
      Kab[[fam]] <- matrix(vcompAB(as.vector(depth_Ref[[fam]]+depth_Alt[[fam]]),as.vector(depth_Ref[[fam]])
                                   ,alpha,epsilon), nrow=nrow(depth_Ref[[fam]]), ncol=ncol(depth_Ref[[fam]]))
    else
      Kab[[fam]] <- bcoef_mat[[fam]]*(1/2)^(depth_Ref[[fam]]+depth_Alt[[fam]])
    Kaa[[fam]] <- bcoef_mat[[fam]]*(1-epsilon)^depth_Ref[[fam]]*epsilon^depth_Alt[[fam]]
    Kbb[[fam]] <- bcoef_mat[[fam]]*(1-epsilon)^depth_Alt[[fam]]*epsilon^depth_Ref[[fam]]
  }
  for(fam in 1:noFam)
    llval = llval + .Call("ll_fs_scaled_err_c",r,Kaa[[fam]],Kab[[fam]],Kbb[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
  return(llval)
}

## r.f.'s are sex-specific
ll_fs_ss_mp_scaled_err <- function(para,depth_Ref,depth_Alt,OPGP,nInd,nSnps,ps,ms,npar,noFam,seqErr,samPara){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- inv.logit2(para[1:npar[1]])
  r[ms,2] <- inv.logit2(para[npar[1]+1:npar[2]])
  if(seqErr)
    epsilon = inv.logit(para[sum(npar)+1])
  else
    epsilon = 0
  if(samPara){
    if(epsilon)
      alpha = exp(para[nSnps+1])
    else
      alpha = exp(para[nSnps])
  }
  
  ## define likelihood
  llval = 0
  # define the density values for the emission probs
  bcoef_mat <- Kab <- Kaa <- Kbb <- vector(mode = "list", length=noFam)
  for(fam in 1:noFam){
    bcoef_mat[[fam]] <- choose(depth_Ref[[fam]]+depth_Alt[[fam]],depth_Ref[[fam]])
    if(samPara)
      Kab[[fam]] <- matrix(vcompAB(as.vector(depth_Ref[[fam]]+depth_Alt[[fam]]),as.vector(depth_Ref[[fam]])
                                   ,alpha,epsilon), nrow=nrow(depth_Ref[[fam]]), ncol=ncol(depth_Ref[[fam]]))
    else
      Kab[[fam]] <- bcoef_mat[[fam]]*(1/2)^(depth_Ref[[fam]]+depth_Alt[[fam]])
    Kaa[[fam]] <- bcoef_mat[[fam]]*(1-epsilon)^depth_Ref[[fam]]*epsilon^depth_Alt[[fam]]
    Kbb[[fam]] <- bcoef_mat[[fam]]*(1-epsilon)^depth_Alt[[fam]]*epsilon^depth_Ref[[fam]]
  }
  for(fam in 1:noFam)
    llval = llval + .Call("ll_fs_ss_scaled_err_c",r,Kaa[[fam]],Kab[[fam]],Kbb[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
  return(llval)
}

## r.f.'s are sex-specific and constrained to the range [0,1] (for unphased data)
ll_fs_up_ss_scaled_err <- function(para,depth_Ref,depth_Alt,bcoef_mat,Kab,config,nInd,nSnps,ps,ms,npar,seqErr){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- inv.logit(para[1:npar[1]])
  r[ms,2] <- inv.logit(para[npar[1]+1:npar[2]])
  if(seqErr)
    epsilon = inv.logit(para[sum(npar)+1])
  else
    epsilon = 0
  ## define likelihood
  llval = 0
  # define the density values for the emission probs
  Kaa <- bcoef_mat*(1-epsilon)^depth_Ref*epsilon^depth_Alt
  Kbb <- bcoef_mat*(1-epsilon)^depth_Alt*epsilon^depth_Ref
  .Call("ll_fs_up_ss_scaled_err_c",r,Kaa,Kab,Kbb,config,nInd,nSnps)
}

compAB <- function(d,a,alpha,ep){
  a1 = rep(0:a,d-a+1)
  a0 = as.vector(sapply(0:(d-a),rep,times=a+1))
  return(sum(exp(lchoose(d-a,a0) + lchoose(a,a1) + (d-a+a1-a0)*log(1-ep) + (a-a1+a0)*log(ep) + 
                   lchoose(d,a) + lbeta(a1+a0+alpha,d-a1-a0+alpha) - lbeta(alpha,alpha))))
}

vcompAB <- Vectorize(compAB, vectorize.args = c("d","a"))

