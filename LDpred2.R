library(data.table)
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(optparse)
library(magrittr)

option_list = list(
	make_option(c("--sst"),action="store",default=FALSE,type="character",help="Input path and filename of summary statistics file"),
	make_option(c("--bfile"),action="store",default=FALSE,type="character",help="Input path and filename of bed file (without .bed)"),
	make_option(c("--out"),action="store",default=FALSE,type="character",help="Input path and filename of output"),
	make_option(c("--ncores"),action="store",default=5,type="integer",help="Input number of cores to run"),
	make_option(c("--geneticPos"),action="store",default=FALSE,type="character",help="Input path for snp_asGeneticPos function's dir argument"))

opt = parse_args(OptionParser(option_list = option_list))
info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))
df = fread(opt$sst,sep=' ')
sumstats <- data.table(df)
sumstats <- sumstats[,c(-5,-(9:11))]
df1 = fread(paste0(opt$bfile,".bim"),sep='\t',header=FALSE)
bim <- data.table(df1)
set1 <- c()
for (i in seq_along(bim$V2)){
	if  (!(bim$V2[i] %in% sumstats$SNP)) {
		set1 <- c(set1,i)
	}
}
bim <- bim[-set1]
sumstats$A0 <- bim$V6
sumstats$SE <- sumstats$SE/sumstats$OR
names(sumstats) <- c("chr","rsid","pos","a1","n_eff","OR","beta_se","p","a0")
sumstats$beta <- log(sumstats$OR)
sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
tmp <- tempfile(tmpdir = ".")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
corr <- NULL
ld <- NULL
fam.order <- NULL
snp_readBed(paste0(opt$bfile,".bed"))
obj.bigSNP <- snp_attach(paste0(opt$bfile,".rds"))
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
info_snp <- snp_match(sumstats, map)
genotype <- obj.bigSNP$genotypes
CHR <- map$chr
POS <- map$pos
POS2 <- snp_asGeneticPos(CHR, POS, dir = opt$geneticPos,ncores=1)
for (chr in 1:1) {
	ind.chr <- which(info_snp$chr == chr)
	ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
	corr0 <- snp_cor(genotype,ind.col=ind.chr2,infos.pos = POS2[ind.chr2],size=3/1000,ncores=opt$ncores)
	if (chr==1){
		ld <- Matrix::colSums(corr0^2)
		corr <- as_SFBM(corr0, tmp)
	}
	else {
		ld <- c(ld, Matrix::colSums(corr0^2))
		corr$add_columns(corr0, nrow(corr))
	}
}
fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order,c("family.ID","sample.ID"),c("FID","IID"))
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(ld,length(ld),chi2=(df_beta$beta / df_beta$beta_se)^2,sample_size=df_beta$n_eff,blocks=NULL)
h2_est <- ldsc[["h2"]]
p_seq <- c(1,0.3,0.1,0.03,0.01,0.003,0.001)
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
grid.param <- expand.grid(p = p_seq,h2 = h2_seq,sparse=FALSE)
beta_grid <- snp_ldpred2_grid(corr, df_beta, grid.param, ncores = opt$ncores)

#if(is.null(obj.bigSNP)){
#	obj.bigSNP <- snp_attach(paste0(opt$bfile,".rds"))
#}
#genotype <- obj.bigSNP$genotypes
#ind.test <- 1:nrow(genotype)
eff <- data.table(info_snp$rsid,info_snp$a1,beta_grid)
write.table(eff,file=opt$out,sep=' ',row.names=FALSE,col.names=FALSE,quote=FALSE)





