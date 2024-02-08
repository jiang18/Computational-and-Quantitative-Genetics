## Software tool
PLINK  
[MPH](https://jiang18.github.io/mph/)

## Data set
The same as in Modules 1 and 2.

## Partioning heritability across chromosomes
```
# The SNP info file is needed by MPH
perl -e 'print "SNP\n"; while(<>){@c=split /\s+/; print "$c[1]\n"}' < geno.qc.bim > snp.info.csv

# Making GRMs
# Input: geno.qc
mkdir chromosomes
for chr in {1..5}
do
    # Extract SNPs from a particular chromosome
    plink --bfile geno.qc --chr $chr --make-bed --out ./chromosomes/$chr
    # Make GRM
    mph --make_grm --binary_genotype ./chromosomes/$chr --snp_info snp.info.csv --num_threads 10 --out ./chromosomes/$chr
    # Remove temporary PLINK files
    rm ./chromosomes/$chr.bim ./chromosomes/$chr.bed ./chromosomes/$chr.fam ./chromosomes/$chr.log ./chromosomes/$chr.nosex
done


# Running REML
# Input: chr.grms.txt and phen.csv
# Make chr.grms.txt
echo "./chromosomes/1 10" > chr.grms.txt
for chr in {2..5}
do
    echo "./chromosomes/$chr 10" >> chr.grms.txt
done

mph --reml --grm_list chr.grms.txt --phenotype phen.csv --trait milk --error_weight milk_wt --num_threads 10 --out ./chromosomes/milk

```
