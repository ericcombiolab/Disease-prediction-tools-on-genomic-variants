require(locfdr)
require(sfsmisc)
library(data.table)
#***************************************************************************************************************************
# For computing empirical Bayes corrected estimates of regression coefficients for polygenic score testing/prediction
#***************************************************************************************************************************

#*******************************************************
# zval: list of z-statistics
# preval: prevalence of disease in the sample
# caseNo: Number of cases
# ctrlNo: Number of controls
# MAF: a vector specifying the minor allele frequencies 
# nulltype, df, bre: parameters for locfdr, please refer to the R package locfdr
#
# Output
# beta.corr.fdrXTweedie: corrected coefficient by tdr*Tweedie method (you can opt to include all variants in the polygenic score; ie choosing of p-value threhsolds can be avoided)
# beta.corr.Tweedie : corrected coefficient by Tweedie method
# beta.corr.fdr : corrected coefficient by tdr weighting method
#*******************************************************

## we assume standardized genotypes
# ********************************************************** 
# R program
# to convert z to corrected Vg with MAF input
# *********************************************************
muaa = 0  ##fixed

D1ss <- function(x, y, xout = x, spar.offset = 0.1384, spl.spar = NULL)
{
  sp <- if (is.null(spl.spar))
  {
    sp <- smooth.spline(x, y)
    smooth.spline(x, y, spar = sp$spar + spar.offset)
  } else smooth.spline(x, y, spar = spl.spar)
  predict(sp, xout, deriv = 1)$y
}

func.RR <- function(RR1, Vg, PA = 0.5)
{
  RR2 = RR1^2
  Paa = (1 - PA)^2
  PAa = 2 * PA * (1 - PA)
  PAA = PA^2
  
  faa = K/(Paa + PAa * RR1 + PAA * RR2)
  fAa = RR1 * faa
  fAA = RR2 * faa
  T = qnorm(1 - faa)  #muaa is set to 0 and residual var set to 1
  muAa = T - qnorm(1 - fAa)
  muAA = T - qnorm(1 - fAA)
  ## overall mean
  mean.all = PAa * muAa + PAA * muAA
  expVg = Paa * (muaa - mean.all)^2 + PAa * (muAa - mean.all)^2 + PAA * 
    (muAA - mean.all)^2
  expVg2 = expVg/(1 + expVg)
  return((expVg2 - Vg)^2)
}

## function to cal. the resulting RR from a given Vg
resRR.func <- function(varexp)
{
  optimize(func.RR, c(1, 1000), Vg = varexp)$minimum
}

## function to cal. power 
power.func <- function(RR1, RR2 = RR1^2, PA = 0.5, K, case.size = caseNo, 
                       ctrl.size = ctrlNo, alpha = 5e-05)
{
  Paa = (1 - PA)^2
  PAa = 2 * PA * (1 - PA)
  PAA = PA^2
  muaa = 0  #fixed
  faa = K/(Paa + PAa * RR1 + PAA * RR2)
  fAa = RR1 * faa
  fAA = RR2 * faa
  
  Paa.case = faa * Paa/K
  PAa.case = fAa * PAa/K
  PAA.case = fAA * PAA/K
  
  Paa.ctrl = (1 - faa) * Paa/(1 - K)
  PAa.ctrl = (1 - fAa) * PAa/(1 - K)
  PAA.ctrl = (1 - fAA) * PAA/(1 - K)
  
  PA.case = PAa.case/2 + PAA.case
  Pa.case = 1 - PA.case
  PA.ctrl = PAa.ctrl/2 + PAA.ctrl
  Pa.ctrl = 1 - PA.ctrl
  
  ORall = PA.case * Pa.ctrl/Pa.case/PA.ctrl

  VarlnOR = 1/(2 * case.size) * (1/PA.case + 1/Pa.case) + 1/(2 * ctrl.size) * 
    (1/PA.ctrl + 1/Pa.ctrl)
  Z = log(ORall)/sqrt(VarlnOR)
  
  
  critR = qnorm(1 - alpha/2)
  critL = qnorm(alpha/2)
  power = 1 - pnorm(critR, mean = Z, sd = 1) + pnorm(critL, mean = Z, 
                                                     sd = 1)
  res.power.func <- list()
  res.power.func$Z = Z
  res.power.func$power = power
  return(res.power.func)
}

muaa = 0  #fixed

## function to be optimized
func.RR <- function(RR1, Vg, PA = 0.5)
{
  RR2 = RR1^2
  Paa = (1 - PA)^2
  PAa = 2 * PA * (1 - PA)
  PAA = PA^2
  
  faa = K/(Paa + PAa * RR1 + PAA * RR2)
  fAa = RR1 * faa
  fAA = RR2 * faa
  T = qnorm(1 - faa)  #muaa is set to 0 and residual var set to 1
  muAa = T - qnorm(1 - fAa)
  muAA = T - qnorm(1 - fAA)
  ## overall mean
  mean.all = PAa * muAa + PAA * muAA
  expVg = Paa * (muaa - mean.all)^2 + PAa * (muAa - mean.all)^2 + PAA * 
    (muAA - mean.all)^2
  expVg2 = expVg/(1 + expVg)
  return((expVg2 - Vg)^2)
}

# ***************************************************************************
# function to cal. Z from a given Vg *
# ***************************************************************************
VgtoZ.func <- function(varexp, PA = 0.5, K = 0.001, case.size = caseNo, 
                       ctrl.size = ctrlNo)
{
  
  RR1 = optimize(func.RR, c(1, 1e+05), Vg = varexp)$minimum
  RR2 = RR1^2
  Paa = (1 - PA)^2
  PAa = 2 * PA * (1 - PA)
  PAA = PA^2
  muaa = 0  #fixed
  faa = K/(Paa + PAa * RR1 + PAA * RR2)
  fAa = RR1 * faa
  fAA = RR2 * faa
  
  Paa.case = faa * Paa/K
  PAa.case = fAa * PAa/K
  PAA.case = fAA * PAA/K
  
  Paa.ctrl = (1 - faa) * Paa/(1 - K)
  PAa.ctrl = (1 - fAa) * PAa/(1 - K)
  PAA.ctrl = (1 - fAA) * PAA/(1 - K)
  
  PA.case = PAa.case/2 + PAA.case
  Pa.case = 1 - PA.case
  PA.ctrl = PAa.ctrl/2 + PAA.ctrl
  Pa.ctrl = 1 - PA.ctrl
  

  Acase = 2 * case.size * PA.case
  acase = 2 * case.size * Pa.case
  Actrl = 2 * ctrl.size * PA.ctrl
  actrl = 2 * ctrl.size * Pa.ctrl
  chisqmat = matrix(c(Acase, acase, Actrl, actrl), nrow = 2, byrow = F)
  Zsq = chisq.test(chisqmat)$statistic
  Z = sqrt(Zsq)
  
  return(Z)
}


power.optim <- function(RR1, obsZ, RR2 = RR1^2, case.size = caseNo, ctrl.size = ctrlNo, 
                        PA = 0.5, K = preval)
{
  Paa = (1 - PA)^2
  PAa = 2 * PA * (1 - PA)
  PAA = PA^2
  muaa = 0  #fixed
  faa = K/(Paa + PAa * RR1 + PAA * RR2)
  fAa = RR1 * faa
  fAA = RR2 * faa
  
  Paa.case = faa * Paa/K
  PAa.case = fAa * PAa/K
  PAA.case = fAA * PAA/K
  
  Paa.ctrl = (1 - faa) * Paa/(1 - K)
  PAa.ctrl = (1 - fAa) * PAa/(1 - K)
  PAA.ctrl = (1 - fAA) * PAA/(1 - K)
  
  PA.case = PAa.case/2 + PAA.case
  Pa.case = 1 - PA.case
  PA.ctrl = PAa.ctrl/2 + PAA.ctrl
  Pa.ctrl = 1 - PA.ctrl
  
  ORall = PA.case * Pa.ctrl/Pa.case/PA.ctrl

  VarlnOR = 1/(2 * case.size) * (1/PA.case + 1/Pa.case) + 1/(2 * ctrl.size) * 
    (1/PA.ctrl + 1/Pa.ctrl)
  Z = log(ORall)/sqrt(VarlnOR)
  
  
  
  return((Z - obsZ)^2)
}

########################################################################################## 
ztoVg.func.bin <- function(obsZ, caseNo, ctrlNo, PA = 0.5, K = 0.1)
{
  
  RR1 = optim(1.13, power.optim, obsZ = obsZ, case.size = caseNo, ctrl.size = ctrlNo, 
              PA = PA, K = K)$par

  RR2 = RR1^2
  
  Paa = (1 - PA)^2
  PAa = 2 * PA * (1 - PA)
  PAA = PA^2
  muaa = 0  #fixed
  
  faa = K/(Paa + PAa * RR1 + PAA * RR2)
  fAa = RR1 * faa
  fAA = RR2 * faa
  
  T = qnorm(1 - faa)  #muaa is set to 0 and residual var set to 1
  muAa = T - qnorm(1 - fAa)
  muAA = T - qnorm(1 - fAA)
  
  ## overall mean
  mean.all = PAa * muAa + PAA * muAA
  Vg = Paa * (muaa - mean.all)^2 + PAa * (muAa - mean.all)^2 + PAA * 
    (muAA - mean.all)^2
  return(Vg/(1 + Vg))
}


# ********************************************************

bias.corr.kernel <- function(zz, xout, bw = "nrd0")
{
  density.obj = density(zz, bw = bw)
  fz.func = splinefun(density.obj$x, density.obj$y)
  Psi.z = log(fz.func(zz)/dnorm(zz))
  truez = D1ss(x = zz, y = Psi.z, xout = xout)
  return(truez)
}

EBay.PRS.binary <- function(zval, preval, caseNo, ctrlNo, MAF, nulltype = 0, 
                            df = 7, bre = 120)
{
  
  # Compute Tweedie formula adjusted z-values and beta
  zall.corr.ker = bias.corr.kernel(zval, xout = zval)
  Vg.corr.ker = mapply(ztoVg.func.bin, obsZ = zall.corr.ker, PA = MAF, 
                       K = preval, caseNo = caseNo, ctrlNo = ctrlNo)
  beta.corr = sqrt(Vg.corr.ker) * sign(zval)
  
  # Compute local fdr
  fdrobj = locfdr(zval, nulltype = nulltype, df = df, bre = bre, plot = 0)
  fdr = fdrobj$fdr
  
  beta.standardized = mapply(ztoVg.func.bin, obsZ = zval, PA = MAF, K = preval, 
                             caseNo = caseNo, ctrlNo = ctrlNo)
  beta.corr.fdr = (1 - fdr) * beta.standardized
  beta.corr.fdrXTweedie = (1 - fdr) * beta.corr
  return(list(beta.corr.fdrXTweedie = beta.corr.fdrXTweedie, beta.corr.Tweedie = beta.corr, 
              beta.corr.fdr = beta.corr.fdr))
  
}

## example
list1=c(0.5)
for (i in 1:10) {
	for (hsq in list1) {
		logistic=read.table(paste('summary statistics file',sep=''),stringsAsFactors=F,header=TRUE,sep=' ')
		OR = logistic[,]
		se = logistic[,]
		selogor = se/OR
		zval = log(OR)/selogor
		logistic$zval <- zval
		MAF.sim = logistic[,]
		obj = EBay.PRS.binary(zval, preval = 0.2, caseNo = num_cases, ctrlNo = num_controls, MAF = MAF.sim, nulltype = 0, df = 7, bre = 120)
		#print(obj)
		logistic$beta <- obj[[2]]
		write.table(logistic,paste('Tweedie_effsize_rep',i,'_',hsq,'.txt',sep=''),row.names=F,col.names=T,quote=F,sep='\t')

}
}
