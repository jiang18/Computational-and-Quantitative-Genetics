## Module 1. Single-variate linear mixed model
https://github.com/jiang18/Computational-and-Quantitative-Genetics/tree/main/module1

### 1. Download data and software
https://www.cog-genomics.org/plink/1.9/  
https://yanglab.westlake.edu.cn/software/gcta/  

### 2. Quality control of genotypes
plink --bfile geno --maf 0.01 --hwe 1e-6 --make-bed --out geno.qc

### 3. Construct GRM
gcta64 --bfile geno.qc --make-grm --threads 4 --out genome

### 4. REML for estimating heritability
gcta64 --reml --grm genome --pheno phen.txt --mpheno 1 --out milk --threads 4

### 5. Alternative REML algorithms
gcta64 --reml --grm genome --pheno phen.txt --mpheno 1 --out milk.fisher --threads 4 --reml-alg 1
gcta64 --reml --grm genome --pheno phen.txt --mpheno 1 --out milk.em --threads 4 --reml-alg 2

### 6. GWAS
gcta64 --mlma --bfile geno.qc --grm genome --pheno phen.txt --mpheno 1 --out milk --threads 4

### 7. LOCO
gcta64 --mlma-loco --bfile geno.qc --pheno phen.txt --mpheno 1 --out milk --threads 14

### 8. Plot
```R
library(qqman)

par(mfrow=c(1,2))

d = read.table("milk.mlma", head=T)
colnames(d)[c(1,3,9)] = c("CHR", "BP", "P")
manhattan(d)

d = read.table("milk.loco.mlma", head=T)
colnames(d)[c(1,3,9)] = c("CHR", "BP", "P")
manhattan(d)

dev.off()

```
