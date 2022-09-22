library(data.table)
library(SuperLearner)

df_valid_score <- data.table()
valid_pheno <- 
df_test_score <- data.table()
test_pheno <- 

sl = SuperLearner(Y=valid_pheno,X=df_valid_score,family=binomial(),SL.library=c("SL.mean","SL.biglasso","SL.ksvm","SL.kernelKnn","SL.earth"))
pred = predict(sl,df_test_score,onlySL=TRUE)
write.table()









