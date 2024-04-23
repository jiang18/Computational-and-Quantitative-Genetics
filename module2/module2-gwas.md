## Software tools
https://github.com/jiang18/slemm  

## Linear regression with phenotypes
```
plink --assoc --bfile geno.qc --allow-no-sex --pheno phen.txt --out milk  
```

```R
library(qqman)

# Manhattan plots and QQ plots
par(mfrow=c(3,2))

d = read.table("milk.mlma", head=T)
colnames(d)[c(1,3,9)] = c("CHR", "BP", "P")
manhattan(d, ylim=c(0,40))
qq(d$P)
lambda = median((d$b/d$se)**2) /  0.455
print(lambda)

d = read.table("milk.loco.mlma", head=T)
colnames(d)[c(1,3,9)] = c("CHR", "BP", "P")
manhattan(d, ylim=c(0,40))
qq(d$P)
lambda = median((d$b/d$se)**2) /  0.455
print(lambda)

d = read.table("milk.qassoc", head=T)
manhattan(d, ylim=c(0,40))
qq(d$P)
lambda = median(d$T**2) /  0.455
print(lambda)

```

## GRAMMAR
```sh
slemm --lmm --phenotype_file phen.csv --bfile geno.qc --trait milk --snp_info_file snp.info.csv --out milk --num_threads 10
perl -e '$_=<>; print "FID $_"; while(<>){print "0 $_"}' < milk.reml.py.txt > milk.residual.txt
plink --assoc --bfile geno.qc --allow-no-sex --pheno milk.residual.txt --pheno-name Py --out milk.residual
```

## GRAMMAR-Gamma
```sh
slemm --lmm --phenotype_file phen.csv --bfile geno.qc --trait milk --snp_info_file snp.info.csv --out milk --num_threads 10
OMP_NUM_THREADS=1 slemm_gamma.py --pfile geno.qc --slemm milk --out milk.gamma.txt
```

## SLEMM-GWA
```sh
slemm --lmm --phenotype_file phen.csv --bfile geno.qc --trait milk --snp_info_file snp.info.csv --out milk --num_threads 10

export OMP_NUM_THREADS=10
for i in `seq 1 5`; do slemm_gwa.py --pfile geno.qc --slemm milk --out milk.chr$i.txt --chr $i; done

mv milk.chr1.txt milk.chrAll.txt
for i in `seq 2 5`; do tail -n +2 milk.chr$i.txt >> milk.chrAll.txt; rm milk.chr$i.txt; done

```
