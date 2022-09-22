library(lassosum)
library(optparse)
library(data.table)
library(penRegSum)
option_list = list(
	make_option(c("--sst"),action="store",default=FALSE,type="character",help="Input path and filename of summary statistics file"),
	make_option(c("--ref"),action="store",default=FALSE,type="character",help="Input path and filename of reference panel"),
	make_option(c("--out"),action="store",default=FALSE,type="character",help="Input path and filename of output"),
	make_option(c("--bfile"),action="store",default=FALSE,type="character",help="Input path and filename of bfile (without .bed)"),
	make_option(c("--lambda"),action="store",default=FALSE,type="double",help="Input lambda"),
	make_option(c("--alpha"),action="store",default=FALSE,type="double",help="Input alpha"))

opt <- parse_args(OptionParser(option_list = option_list))
df1 <- fread(opt$sst,sep=' ')
ss <-data.table(df1)
ss <- ss[!P == 0]
cor <- p2cor(p = ss$P,n = ss$NMISS,sign=log(ss$OR))

elast <- elastSum(cor,opt$ref,opt$lambda,opt$alpha)

eff <- data.table(ss$SNP,ss$A1,elast$beta)
write.table(eff,file=opt$out,sep=' ',row.names=FALSE,col.names=FALSE,quote=FALSE)
