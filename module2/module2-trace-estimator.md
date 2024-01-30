## The demo is based on the dataset that is also used in Module 1.

## Construct the genomic relationship matrix using MPH (https://jiang18.github.io/mph/)
```
perl -e 'print "SNP\n"; while(<>){@c=split /\s+/; print "$c[1]\n"}' < geno.qc.bim > snp.info.csv
mph --make_grm --binary_genotype geno.qc --snp_info snp.info.csv --output demo2
```

## R script for comparing exact trace with stochastic estimates
```R
source("mph_functs.R")

grm = read_grm("demo2")
# Note the upper triangle is empty and needs to be filled.
grm = grm + t(grm)
diag(grm) = diag(grm)/2

exact_tr = sum(diag(grm))

# sample size
n = nrow(grm)
# number of replicates
nrep = 20

est_tr = matrix(nrow=nrep, ncol=5)
colnames(est_tr) = seq(100,500,100)

# The number of random vectors is set to 100, 200, 300, 400, and 500.
for(i in 1:5) {
    nrand = i*100
    for(j in 1:nrep) {
        cur_tr = 0
        for(k in 1:nrand) {
            x = rnorm(n=n)
            cur_tr = cur_tr + sum(x * (grm %*% x))
        }
        cur_tr = cur_tr/nrand
        
        est_tr[j, i] = cur_tr
        print(paste(j, "replicates completed."))
    }
}

```