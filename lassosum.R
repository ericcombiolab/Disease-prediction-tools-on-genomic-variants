library(optparse)
library(lassosum)
library(data.table)
library(methods)
library(magrittr)
library(parallel)

option_list = list(
	 make_option(c("--sst"),action="store",default=FALSE,type="character",help="Input path and filename of summary statistics file"),
	 make_option(c("--bfile"),action="store",default=FALSE,type="character",help="Input path and filename of bed file (without .bed)"),
	 make_option(c("--out"),action="store",default=FALSE,type="character",help="Input path and filename of output"),
	 make_option(c("--ncores"),action="store",default=5,type="integer",help="Input number of cores to run"),
	 make_option(c("--covariate"),action="store",default=FALSE,type="character",help="Input path and filename of covariate"))

opt = parse_args(OptionParser(option_list = option_list))
cl <- makeCluster(opt$ncores)
df1 <- fread(opt$sst,sep=' ')
covariate <- fread(opt$covariate)
cov <- as.data.frame(covariate)
ss <-data.table(df1)
ss <- ss[!P == 0]
target.pheno <- fread(paste0(opt$bfile,".fam"),sep=' ')[,c(1,2,6)]
names(target.pheno) <- c("FID","IID","PHENO")
ld.file <- "EUR.hg19"
cor <- p2cor(p = ss$P,n = ss$NMISS,sign=log(ss$OR))
fam <- fread(paste0(opt$bfile, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]
out <- lassosum.pipeline(cor = cor,chr=ss$CHR,pos = ss$BP,snp=ss$SNP,A1 = ss$A1,A2 = ss$A2,ref.bfile = opt$bfile,test.bfile = opt$bfile,sample=5000,LDblocks = ld.file, cluster=cl)
v <- validate(out,pheno=as.data.frame(target.pheno),covar=cov)
eff <- data.table(out$sumstats$snp,out$sumstats$A1,v$best.beta)
write.table(eff,file=opt$out,sep=' ',row.names=FALSE,col.names=FALSE,quote=FALSE)




