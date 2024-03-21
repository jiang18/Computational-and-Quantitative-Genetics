# Data  
- [Wheat](https://cran.r-project.org/web/packages/BLR/index.html)
- 599 individuals and 1279 markers
- Four traits
  - T1 (grain yield in environment 1)
  - T2 (grain yield in environment 2)
  - T3 (grain yield in environment 3)
  - T4 (grain yield in environment 4)
- PLINK files are provided in the current directory.
  - Note that the SNP and pedigree information is not filled in bim/fam. 

# Split the data into training/validation
```sh
infile="wheat/wheat.csv"
shuffled="shuffled.csv"
(tail -n +2 $infile | shuf) > $shuffled

(head -1 $infile && head -500 $shuffled) > training.csv
(head -1 $infile && tail -99 $shuffled) > validation.csv

```

# Run SLEMM
```sh
# generate the snp info file.
perl -e 'print "SNP\n"; while(<>){@c=split /\s+/; print "$c[1]\n"}' < wheat/wheat.bim > snp.info.csv

# run REML with SLEMM
# --iter_weighting specifies iterative SNP weighting.
trt=T1
slemm --reml --iter_weighting --bfile wheat/wheat --snp_info snp.info.csv --phenotype training.csv --trait $trt --out $trt
# compute GEBVs
slemm --pred --snp_estimate $trt.reml.snp.csv --bfile wheat/wheat --output $trt.gebv.csv

```

# Run LDAK
```sh
# prepare phenotype files
trt=T1
perl -e '$_=<>; while(<>) {@c=split /,/; print "0 $c[0] $c[1]\n";}' < training.csv > $trt.training.pheno
perl -e '$_=<>; while(<>) {@c=split /,/; print "0 $c[0] $c[1]\n";}' < validation.csv > $trt.validation.pheno

# prepare heritability option file
perl -e 'while(<>){@c=split /\s+/; print "$c[1] 2e-4\n"}' < wheat/wheat.bim > opt.hers

# run BOLT
ldak --bolt bolt.$trt.opt --bfile wheat/wheat --ind-hers opt.hers --pheno $trt.training.pheno --cv-proportion .1
ldak --calc-scores bolt.$trt.opt --bfile wheat/wheat --scorefile bolt.$trt.opt.effects --power 0 --pheno $trt.validation.pheno

# run BayesR
ldak --bayesr bayesr.$trt.opt --bfile wheat/wheat --ind-hers opt.hers --pheno $trt.training.pheno --cv-proportion .1
ldak --calc-scores bayesr.$trt.opt --bfile wheat/wheat --scorefile bayesr.$trt.opt.effects --power 0 --pheno $trt.validation.pheno

```

# Run MCMC-BayesR
```sh
# prepre input files

cp -r wheat mcmc-bayesr
```

```R
val = read.csv("validation.csv")
pheno = read.csv("mcmc-bayesr/wheat.csv")
pheno[pheno$ID %in% val$ID, -1] = NA

fam = read.table("mcmc-bayesr/wheat.fam")
fam = fam[,1:5]
fam = cbind(fam, pheno[,-1])

write.table(fam, file="mcmc-bayesr/wheat.fam", quote=F, na="NA", row.names=F, col.names=F)

```

```sh
bayesRv2 -bfile mcmc-bayesr/wheat -n 1 -blocksize 4 -msize 500 -permute -nthreads 4 -out T1 -numit 50000 -burnin 20000 -vara 0.001 -vare 0.55 
bayesRv2 -bfile mcmc-bayesr/wheat -n 1 -out T1.pred -predict -model T1.model -freq T1.frq -param T1.param

```



# Compute metrics
```R
val = read.csv("validation.csv")

trt = "T1"
pred = read.csv(paste0(trt, ".gebv.csv"))
colnames(pred)[1] = "ID"

merged = merge(pred, val, by="ID")

cor(merged$GEBV, merged[[trt]])

```
