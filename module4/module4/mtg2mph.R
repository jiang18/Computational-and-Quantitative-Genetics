source("../demo/mph_functs.R")
library(data.table)

fam = fread("tst.fam")
mtg = fread("tst.grm")
ss = nrow(fam)

iid = fam[[2]]
grm = matrix(nrow=ss, ncol=ss)
grm[upper.tri(grm, diag = TRUE)] <- mtg[[3]]
grm = t(grm)

write_grm("example1-3-2", iid, grm)

######
dat = fread("tst.dat")
dat = dat[,-1]
colnames(dat) = c("IID", 1:2)
fwrite(dat, file="example1-3-2.pheno.csv", sep=",", quote=F)
