library(rMVP)
library(simer)
library(data.table)
library(optparse)

option_list = list(
        make_option(c("--num"),action="store",default=FALSE,type="integer",help="Input replication number")
)

opt = parse_args(OptionParser(option_list = option_list))
popu = 'EUR'
chro = 22
sim_num = opt$num
MVP.Data(fileBed=paste('./',popu,'/hapmap3_',popu,'_chr',chro,sep=''),filePhe=NULL,fileKin=FALSE,filePC=FALSE,out=paste('./',popu,'/mvp.plink_',sim_num,sep=''))

num_snps = 20

print('Loading genotypes')
pop.geno = bigmemory::attach.big.matrix(paste('./',popu,'/mvp.plink_',sim_num,'.geno.desc',sep=''))
pop.map = read.table(paste('./',popu,'/mvp.plink_',sim_num,'.geno.map',sep=''),header=TRUE)

SP <- param.annot(pop.map = pop.map, qtn.num = list(tr1 = num_snps), qtn.model = "A + D + A:D") 
SP <- param.geno(SP = SP, pop.geno = pop.geno)

#SP <- param.geno(SP = SP, pop.marker = 2e4, pop.ind = 1e2)

SP <- param.pheno(
    SP = SP,
    pop.ind = 90000,
    phe.model = list(tr1 = "T1 = A + D + A:D + E"),
    phe.h2A = list(tr1 = 0.3),
    phe.h2D = list(tr1 = 0.1),
    phe.h2GxG = list(tr1 = list("A:D" = 0.1))
)
print('Annotation')
SP <- annotation(SP)
print('genotype')
SP <- genotype(SP)
print('phenotype')
SP <- phenotype(SP)    

SP <- param.simer(
  SP = SP,
  out = paste(popu,'_chr',chro,'_rep',sim_num,sep=''),
  replication = sim_num,
  seed.sim = 1,
  outpath = paste('./',popu,'/',sep=''),
  out.format = "plink",
  ncpus = 30
)


# Run Simer
SP <- simer(SP)

