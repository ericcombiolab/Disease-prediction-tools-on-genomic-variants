library(EBPRS)
library(optparse)
library(data.table)

option_list = list(
	make_option(c("--sst"),action="store",default=FALSE,type="character",help="Input path and filename of summary statistics file"),
	make_option(c("--out"),action="store",default=FALSE,type="character",help="Input path and filename of output"),
	make_option(c("--num_cases"),action="store",default=FALSE,type="integer",help="Input the number of cases"),
	make_option(c("--num_controls"),action="store",default=FALSE,type="integer",help="Input the number of controls"))

opt = parse_args(OptionParser(option_list = option_list))
df1 <- fread(opt$sst,sep=' ')
ss <-data.table(df1)
ss <- ss[,c(-1,-3,-6,-8,-9,-(11:14))]
setcolorder(ss,c("A1","A2","OR","P","SNP"))

out <- EBPRS(train=ss,N1=opt$num_cases,N0=opt$num_controls)
eff <- data.table(out$result$SNP,out$result$A1,out$result$effectsize)
write.table(eff,file=opt$out,sep=' ',row.names=FALSE,col.names=FALSE,quote=FALSE)
