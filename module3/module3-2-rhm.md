## Software tools
PLINK
[MPH](https://jiang18.github.io/mph/)

## Data set
The same as in Modules 1 and 2.

## Regional heritability mapping
Sliding windows: size = 5 Mbp and moving = 2.5 Mbp

```
# Making a GRM from whole-genome SNPs
mph --make_grm --binary_genotype geno.qc --snp_info snp.info.csv --output genome

# RHM
mkdir rhm
# GRM lists
echo "genome" > deduct.list.txt
echo "./rhm/region" >> deduct.list.txt

echo "./rhm/remaining" > remaining.list.txt

echo "./rhm/region" > rhm.list.txt
echo "./rhm/remaining"  >> rhm.list.txt

for chr in {1..5}; do
    for wi in {0..39}; do
        # Extract SNPs from a particular region 
        wstart=$(echo "$wi * 2.5" | bc)
        wend=$(echo "$wstart + 5" | bc)
        plink --bfile geno.qc --chr $chr --from-mb $wstart --to-mb $wend --make-bed --out ./rhm/region
        # Make a GRM for a particular region
        mph --make_grm --binary_genotype ./rhm/region --snp_info snp.info.csv --num_threads 10 --out ./rhm/region
        # Make a GRM from the remaining SNPs
        mph --deduct_grm --grm_list deduct.list.txt --out ./rhm/remaining
        # REML for the full model and the reduced model
        mph --reml --grm_list rhm.list.txt --phenotype phen.csv --trait milk --num_threads 10 --out ./rhm/$chr.$wi.full
        mph --reml --grm_list remaining.list.txt --phenotype phen.csv --trait milk --num_threads 10 --out ./rhm/$chr.$wi.reduced
    done
done

```

## Summarize the REML logLL
```R
lrt = data.frame(SNP=c(1:200), CHR=rep(NA, 5*40), BP=rep(NA, 5*40), CHISQ=rep(NA, 5*40), P=rep(NA, 5*40))
k = 0
for (chr in 1:5) {
    for (i in 0:39) {
        k = k + 1

        full = read.csv( paste0("./rhm/", chr, ".", i, ".full.mq.iter.csv") )
        reduced = read.csv( paste0("./rhm/", chr, ".", i, ".reduced.mq.iter.csv") )

        lrt[k, 2] = chr
        lrt[k, 3] = (i+1) * 2.5e6
        lrt[k, 4] = 2*(full$logLL[nrow(full)] - reduced$logLL[nrow(reduced)])
        lrt[k, 5] = pchisq(lrt[k, 4], df=1, lower=F) / 2
    }
}

library(qqman)

par(mfrow=c(2,1))
manhattan(lrt, genomewide=-log10(0.05/nrow(lrt)), suggestive=F, main="Regional heritability mapping")

d = read.table("milk.mlma", head=T)
colnames(d)[c(1,3,9)] = c("CHR", "BP", "P")
manhattan(d, genomewide=-log10(0.05/nrow(d)), suggestive=F, main="Single-marker mixed-model associations")

```
